# -*- coding: utf-8 -*-
'''网页界面: 上传 STEP（或用内置结构生成）、核对面类型和多孔体、填工况，提交到后台排队自动计算（第 2–8 步）。

本机使用（在仓库根目录）:
    pip install gmsh streamlit plotly
    streamlit run app.py
多人使用（服务器）: 见 docs/11-多人网页部署.md，用 Dockerfile / docker-compose.yml 部署。

计算不在网页进程里做: 点“提交计算”后任务进入队列（job_queue.py），由后台进程 worker.py 逐个计算，
网页关掉也不影响；网页发现后台进程没在运行时会自动启动它。每个用户只看得到自己的任务。

环境变量（都可以不设，不设就是单人本机模式）:
    ZXDCFD_WORKDIR              工作根，默认 仓库/运行目录
    ZXDCFD_USERS                用户和访问码，格式 “张三:码1,李四:码2”；设了就要登录
    ZXDCFD_USERS_FILE           或者用 JSON 文件 {"张三": "码1"}
    ZXDCFD_ADMINS               管理员用户名，逗号分隔，能看到所有人的任务
    ZXDCFD_RUN_MODE             固定运行方式（容器 / 本机），设了用户就不能选；服务器镜像里是 本机
    ZXDCFD_MAX_CELLS            单元数上限的最大值，默认 400000
    ZXDCFD_MAX_TIME             单次求解时限的最大值（秒），默认 7200
    ZXDCFD_MAX_JOBS_PER_USER    每人同时排队 + 计算的任务数，默认 2
    ZXDCFD_START_WORKER         设为 0 时网页不自动启动后台进程（由别的方式启动 worker.py）
自动返工和网格无关性检查由 auto_run.py 完成，规则见 docs/09-网页界面与自动运行.md。
'''
import hashlib
import hmac
import json
import os
import shutil
import subprocess
import threading
import time
import uuid

import gmsh
import pandas as pd
import streamlit as st

import auto_run
import job_queue
import ui_tools
import worker
from generate_boundary import 识别几何
from generate_case import 默认步数
from generate_structure import generate_step, 默认参数, 结构清单, 结构标注
from run_case import 默认时限, 默认镜像

# Streamlit 在子线程里运行脚本，gmsh.initialize() 默认要装信号处理，只能在主线程，所以这里先用
# interruptible=False 初始化；之后第 1–4 步的函数发现已初始化就直接复用，不再初始化和 finalize。
if not gmsh.isInitialized():
    gmsh.initialize(interruptible=False)

来源选项 = ('上传 STEP 文件', '用内置结构生成')
页面选项 = ('新建计算', '我的任务')
状态图标 = {'排队': '排队中', '运行中': '计算中', '完成': '已完成', '出错': '出错', '已取消': '已取消'}

st.set_page_config(page_title='zxdcfd 三维换热计算', layout='wide')
状态 = st.session_state


@st.cache_resource
def _gmsh锁对象():
    return threading.Lock()


def gmsh锁():
    '''gmsh 是进程内全局状态，多个用户的会话在同一进程的不同线程里，读 STEP、生成几何、三角化要排队做。'''
    return _gmsh锁对象()


def 读用户表():
    表 = {}
    文件 = os.environ.get('ZXDCFD_USERS_FILE')
    if 文件:
        with open(文件, encoding='utf-8') as f:
            表.update({str(k): str(v) for k, v in json.load(f).items()})
    for 项 in (os.environ.get('ZXDCFD_USERS') or '').split(','):
        if ':' in 项:
            名, 码 = 项.split(':', 1)
            表[名.strip()] = 码.strip()
    return {k: v for k, v in 表.items() if k and v}


用户表 = 读用户表()
管理员 = {n.strip() for n in (os.environ.get('ZXDCFD_ADMINS') or '').split(',') if n.strip()}
固定运行方式 = os.environ.get('ZXDCFD_RUN_MODE') or None
限额 = {'单元': int(os.environ.get('ZXDCFD_MAX_CELLS') or auto_run.默认单元上限),
        '时限': int(os.environ.get('ZXDCFD_MAX_TIME') or 7200),
        '每人任务': int(os.environ.get('ZXDCFD_MAX_JOBS_PER_USER') or 2)}


def 新工作目录():
    目录 = os.path.join(job_queue.工作根(), '会话', time.strftime('%Y%m%d-%H%M%S') + '-' + uuid.uuid4().hex[:6])
    os.makedirs(目录, exist_ok=True)
    return 目录


def 换几何(step路径, 标注=None, 说明=''):
    '''读入新的 STEP: 自动识别，套用标注。'''
    with gmsh锁():
        几何 = 识别几何(step路径)
    面类 = {f['编号']: f['类型'] for f in 几何['面']}
    if 标注:
        for t, 类 in 标注.get('面', {}).items():
            面类[int(t)] = 类
    状态.update(step路径=step路径, 几何=几何, 自动面类={f['编号']: f['类型'] for f in 几何['面']},
              面类=面类, 多孔体=list(标注['多孔体']) if 标注 and '多孔体' in 标注 else list(几何['多孔体']),
              几何说明=说明, 边界版本=状态.get('边界版本', 0) + 1)
    状态.pop('三角网格', None)


def 指纹(*对象):
    return hashlib.sha1(json.dumps(对象, ensure_ascii=False, sort_keys=True, default=str).encode()).hexdigest()


def 可以看(任务):
    return 任务.get('用户') == 状态.用户 or 状态.用户 in 管理员


def 时间(t):
    return time.strftime('%m-%d %H:%M', time.localtime(t)) if t else ''


def 用时(任务):
    if not 任务.get('开始时间'):
        return ''
    秒 = (任务.get('结束时间') or time.time()) - 任务['开始时间']
    return '%d 分 %02d 秒' % (秒 // 60, 秒 % 60)


# ================================================================ 登录
st.title('zxdcfd 三维换热计算')
if not 用户表:
    状态.用户 = '本机'
elif not 状态.get('用户'):
    st.caption('请用管理员给你的用户名和访问码登录。')
    with st.form('登录'):
        名 = st.text_input('用户名', key='登录名')
        码 = st.text_input('访问码', type='password', key='登录码')
        登录 = st.form_submit_button('登录', type='primary')
    if 登录:
        if 名 in 用户表 and hmac.compare_digest(码.encode(), 用户表[名].encode()):
            状态.用户 = 名
            st.rerun()
        time.sleep(1)                     # 放慢猜码
        st.error('用户名或访问码不对')
    st.stop()

if os.environ.get('ZXDCFD_START_WORKER') != '0':
    worker.确保运行()

# 切到“我的任务”页时新建页的控件不显示，Streamlit 会丢掉它们的值；每次运行先把值写回，切回来时保持原样。
for 键 in list(状态.keys()):
    if 键 in ('几何来源', '结构号', '网格尺寸', '步数', '单元上限', '最多返工', '网格检查', '对照一维', '运行方式',
              '时限', '任务说明') or 键.startswith(('参数_', '工况_')):
        状态[键] = 状态[键]

# ================================================================ 侧栏
if 状态.pop('跳转任务', None):
    状态.页面 = '我的任务'
参数任务 = st.query_params.get('任务')
if 参数任务 and 状态.get('_已读参数') != 参数任务:
    状态._已读参数 = 参数任务
    状态.页面, 状态.选中任务 = '我的任务', 参数任务
with st.sidebar:
    st.markdown('**%s**' % ('单机模式' if not 用户表 else '用户：' + 状态.用户))
    if 用户表 and st.button('退出登录', key='退出'):
        for 键 in list(状态.keys()):
            del 状态[键]
        st.rerun()
    页面 = st.radio('页面', 页面选项, key='页面')
    全部 = job_queue.列出()
    st.caption('计算服务：%s；排队 %d 个，计算中 %d 个'
               % ('运行中' if worker.在运行() else '未启动',
                  sum(t['状态'] == '排队' for t in 全部), sum(t['状态'] == '运行中' for t in 全部)))


# ================================================================ 结果显示（我的任务页用）
def 显示结果(总结, 号):
    显示 = {'通过': st.success, '网格无关性未确认': st.warning, '需要人工': st.error, '未通过': st.error}
    显示.get(总结['结论'], st.info)('**%s**　%s' % (总结['结论'], 总结['说明']))

    最终 = 总结.get('最终')
    if 最终 and 最终.get('结果'):
        结果 = 最终['结果']
        c = st.columns(5)
        c[0].metric('进出口压差', '%.4g Pa' % 结果['进出口压差'])
        c[1].metric('出口温度（流量加权）', '%.4g K' % 结果['出口温度'])
        c[2].metric('吸热量', '%.4g MW' % (结果['吸热量'] / 1e6))
        c[3].metric('入口流量', '%.4g m³/s' % 结果['入口流量'])
        c[4].metric('网格', '%d 个单元' % 最终['单元数'], '尺寸 %.3g m' % 最终['网格尺寸'], delta_color='off')
        if 结果.get('多孔压降估算') is not None:
            st.caption('一维估算（均匀流过多孔区，仅供参考）: 多孔压降 %.4g Pa，出口温度 %.4g K。'
                       '进出口压差是静压差，入口和出口流速不同时含动压变化，不能直接和一维压降比。'
                       % (结果['多孔压降估算'], 结果['出口温度估算']))
        st.subheader('第 7 步检查')
        st.dataframe([{'项目': 项['项目'], '结果': '通过' if 项['通过'] else '未通过', '说明': 项['说明']}
                      for 项 in 结果['检查']], hide_index=True, width='stretch')

    if 总结['网格无关性']:
        st.subheader('网格无关性')
        st.dataframe([{'第几次': 行['编号'], '网格尺寸 m': round(行['网格尺寸'], 4), '单元数': 行['单元数'],
                       '进出口压差 Pa': round(行['进出口压差'], 4), '出口温度 K': round(行['出口温度'], 3),
                       '吸热量 MW': round(行['吸热量'] / 1e6, 5),
                       '压差变化': '%.2f%%' % (100 * 行['相对变化']['进出口压差']) if '相对变化' in 行 else '',
                       '吸热量变化': '%.2f%%' % (100 * 行['相对变化']['吸热量']) if '相对变化' in 行 else ''}
                      for 行 in 总结['网格无关性']], hide_index=True, width='stretch')

    st.subheader('每一次运行')
    st.dataframe([{'第几次': r['编号'], '网格层': r['层'], '网格尺寸 m': round(r['网格尺寸'], 4) if r['网格尺寸'] else None,
                   '单元数': r['单元数'], '步数': r['步数'], '完成步数': r['完成步数'],
                   '第 8 步': r['结论']['结论'], '自动修改 / 原因': r['改法'] or r['异常'] or ''}
                  for r in 总结['尝试']], hide_index=True, width='stretch')

    if 最终 and 最终.get('算例目录'):
        历史 = ui_tools.读历史(最终['算例目录'])
        try:
            import plotly.graph_objects as go
            c = st.columns(2)
            if 历史['残差']:
                图 = go.Figure([go.Scatter(y=值, name=名, mode='lines') for 名, 值 in 历史['残差'].items()])
                图.update_layout(title='残差（每步初始残差）', yaxis_type='log', height=320,
                                 margin={'l': 10, 'r': 10, 't': 40, 'b': 10}, xaxis_title='步')
                c[0].plotly_chart(图, width='stretch')
            if 历史['监测'].get('outletT'):
                图 = go.Figure([go.Scatter(y=历史['监测']['outletT'], name='出口温度 K', mode='lines')])
                if 历史['监测'].get('inletP'):
                    图.add_trace(go.Scatter(y=历史['监测']['inletP'], name='入口 p/ρ (m²/s²)', mode='lines',
                                            yaxis='y2'))
                图.update_layout(title='监测量', height=320, margin={'l': 10, 'r': 10, 't': 40, 'b': 10},
                                 xaxis_title='步', yaxis={'title': 'K'},
                                 yaxis2={'title': 'm²/s²', 'overlaying': 'y', 'side': 'right'})
                c[1].plotly_chart(图, width='stretch')
        except ImportError:
            pass

    c = st.columns(3)
    if 最终 and 最终.get('算例目录') and os.path.isdir(最终['算例目录']):
        包键 = '算例包_' + 号
        if 包键 in 状态:
            c[0].download_button('下载最终算例（zip）', 状态[包键], file_name='算例-%s.zip' % 号,
                                 mime='application/zip', key='下载算例')
        elif c[0].button('打包最终算例', key='打包算例'):
            状态[包键] = ui_tools.打包(最终['算例目录'])
            st.rerun()
    路径 = os.path.join(总结['输出目录'], auto_run.报告文件名)
    if os.path.isfile(路径):
        with open(路径, encoding='utf-8') as f:
            c[1].download_button('下载运行报告', f.read(), file_name='报告-%s.txt' % 号, key='下载报告')
    c[2].download_button('下载总结（json）', ui_tools.转json(总结), file_name=auto_run.总结文件名,
                         mime='application/json', key='下载总结')

    with st.expander('这个结果能信到什么程度', expanded=总结['结论'] == '通过'):
        st.markdown(
            '- **自动保证的**：求解跑完、残差和监测值稳定、进出口流量守恒、入口流量与设定一致；'
            '勾选了网格无关性时，再加密一层网格结果变化在容差内（压差 3%、吸热量 2%）。\n'
            '- **不能自动保证的**：物理模型本身。现在是不可压、密度取入口温度下的值，层流加常数湍流粘度，'
            '外壁绝热，多孔阻力和换热系数由你填写。转弯、分离、导流板这类问题对湍流模型和近壁处理敏感，'
            '这个模型只能看趋势，不能当设计值。\n'
            '- **结果准不准，最终要和实测或可靠的设计计算对比**；数值检查通过只说明“按这个模型算对了”。')


def 任务面板(号):
    任务 = job_queue.读(号)
    c = st.columns([2, 1, 1, 1])
    c[0].markdown('**任务 %s**　%s' % (号, 任务.get('说明') or ''))
    文字 = 状态图标[任务['状态']]
    if 任务['状态'] == '排队':
        位置 = job_queue.排队位置(号)
        文字 += '（前面还有 %d 个）' % 位置 if 位置 else '（下一个就轮到）'
    c[1].markdown('状态：**%s**' % 文字)
    c[2].markdown('用时：%s' % (用时(任务) or '-'))
    if 任务['状态'] in job_queue.未结束:
        if c[3].button('取消任务', key='取消_' + 号):
            job_queue.取消(号)
            st.rerun()
    elif c[3].button('删除任务', key='删除_' + 号):
        job_queue.删除(号)
        状态.选中任务 = None
        st.rerun()
    if 任务['状态'] in ('出错', '已取消'):
        st.error('%s：%s' % (任务['状态'], 任务.get('结论说明') or ''))
    日志 = job_queue.日志(号, 末尾行数=400)
    with st.expander('计算过程', expanded=任务['状态'] in job_queue.未结束 or 任务['状态'] == '出错'):
        st.code(日志 or '（还没有开始）', language=None)
    return 任务


@st.fragment(run_every=3)
def 进行中面板(号):
    任务 = 任务面板(号)
    if 任务['状态'] not in job_queue.未结束:
        st.rerun(scope='app')


if 页面 == '我的任务':
    st.header('我的任务' if 状态.用户 not in 管理员 else '全部任务（管理员）')
    列表 = [t for t in 全部 if 可以看(t)]
    if not 列表:
        st.info('还没有任务。到“新建计算”提交一个。')
        st.stop()
    表 = [{'任务号': t['号'], **({'用户': t['用户']} if 状态.用户 in 管理员 else {}), '说明': t.get('说明') or '',
          '状态': 状态图标[t['状态']], '结论': t.get('结论') or '', '提交': 时间(t['提交时间']), '用时': 用时(t)}
         for t in 列表]
    st.dataframe(表, hide_index=True, width='stretch')
    号表 = [t['号'] for t in 列表]
    if 状态.get('选中任务') in 号表:
        状态['选任务'] = 状态.选中任务
    号 = st.selectbox('查看任务', 号表, key='选任务',
                     format_func=lambda h: '%s %s %s' % (h, 状态图标[next(t for t in 列表 if t['号'] == h)['状态']],
                                                         next(t for t in 列表 if t['号'] == h).get('说明') or ''))
    状态.选中任务 = 号
    st.query_params['任务'] = 号
    状态._已读参数 = 号
    if 参数任务 and 参数任务 not in 号表:
        st.warning('链接里的任务 %s 不存在或不属于你。' % 参数任务)
    任务 = job_queue.读(号)
    if 任务['状态'] in job_queue.未结束:
        进行中面板(号)
    else:
        任务面板(号)
        总结 = job_queue.总结(号)
        if 任务['状态'] == '完成' and 总结:
            st.header('结果')
            显示结果(总结, 号)
    st.caption('任务结束 %s 天后自动删除；可以把本页链接收藏，之后登录回来看。'
               % (os.environ.get('ZXDCFD_KEEP_DAYS') or 7))
    st.stop()

# ================================================================ 新建计算
if '任务' in st.query_params:
    del st.query_params['任务']
st.caption('上传几何 → 核对边界 → 填工况 → 提交。后台自动计算，出错时按规则修改网格、步数或时限后重跑，'
           '通过后再加密网格确认结果不随网格变化。')
if '工作目录' not in 状态:
    状态.工作目录 = 新工作目录()
if '工况' not in 状态:
    状态.工况, 状态.工况版本 = ui_tools.默认工况(), 0

# ================================================================ 1 几何
st.header('1  几何')
来源 = st.radio('几何来源', 来源选项, horizontal=True, key='几何来源')
if 来源 == 来源选项[0]:
    文件 = st.file_uploader('STEP 文件（.step / .stp，单位按 m 处理；多个实体会自动 fragment，交界面共享）',
                            type=['step', 'stp'], key='上传STEP')
    if 文件 is not None:
        标识 = 指纹(文件.name, 文件.size)
        if 状态.get('上传标识') != 标识:
            路径 = os.path.join(状态.工作目录, '上传', 'model.step')
            os.makedirs(os.path.dirname(路径), exist_ok=True)
            with open(路径, 'wb') as f:
                f.write(文件.getbuffer())
            try:
                换几何(路径, 说明='上传: %s' % 文件.name)
                状态.上传标识 = 标识
            except Exception as 错:
                st.error('STEP 读不了: %s' % 错)
else:
    列 = st.columns([1, 3])
    结构号 = 列[0].selectbox('结构', list(结构清单), key='结构号',
                            format_func=lambda k: '%s %s' % (k, 结构清单[k]['名称']))
    参数 = {}
    参数列 = 列[1].columns(4)
    for i, (名, 值) in enumerate(默认参数(结构号).items()):
        参数[名] = 参数列[i % 4].number_input('%s (m)' % 名, value=float(值), min_value=0.001, format='%.4g',
                                             key='参数_%s_%s' % (结构号, 名))
    if st.button('生成几何', key='生成几何'):
        try:
            目录 = os.path.join(状态.工作目录, '生成')
            with gmsh锁():
                路径 = generate_step(结构号, 参数, 目录, 'model')
                标注 = 结构标注(结构号, 参数, 路径)
            换几何(路径, 标注, 说明='内置结构 %s %s' % (结构号, 结构清单[结构号]['名称']))
        except Exception as 错:
            st.error('生成失败: %s' % 错)

if not 状态.get('step路径'):
    st.info('先上传 STEP 文件，或选一种内置结构生成几何。')
    st.stop()

几何 = 状态.几何
st.success('%s。体 %d 个，面 %d 个，包围盒对角线 %.4g m。'
           % (状态.几何说明, len(几何['体']), len(几何['面']), 几何['对角线']))

# ================================================================ 2 边界
st.header('2  边界')
st.caption('面类型已按规则自动猜好（两个体共用的面是 interface；模型最高的水平面是 inlet，最低的是 outlet；'
           'x=0 且是最小 x 的平面是 symmetry；其余 walls）。请对照三维图逐个核对，猜错的直接在表里改。'
           '入口在侧面、转弯入口等情况一定要改。')
外表面 = [f for f in sorted(几何['面'], key=lambda f: f['编号']) if 状态.自动面类[f['编号']] != 'interface']
内部面 = [f['编号'] for f in 几何['面'] if 状态.自动面类[f['编号']] == 'interface']

左, 右 = st.columns([3, 2])
with 右:
    表 = [{'面': f['编号'], '类型': 状态.面类[f['编号']], '自动识别': 状态.自动面类[f['编号']],
          '面积 m²': round(f['面积'], 6), 'x': round(f['中心'][0], 4), 'y': round(f['中心'][1], 4),
          'z': round(f['中心'][2], 4)} for f in 外表面]
    改后 = st.data_editor(
        pd.DataFrame(表), key='面表_%d' % 状态.边界版本, hide_index=True, width='stretch',
        disabled=['面', '自动识别', '面积 m²', 'x', 'y', 'z'],
        column_config={'类型': st.column_config.SelectboxColumn('类型', options=list(ui_tools.可选面类型),
                                                               required=True)})
    for 行 in 改后.to_dict('records'):
        状态.面类[int(行['面'])] = 行['类型']
    if 内部面:
        st.caption('内部面（interface，不能设为边界）: %s' % ', '.join(map(str, sorted(内部面))))
    体标签 = {v['编号']: '体 %d（体积 %.4g m³，中心 z=%.4g）' % (v['编号'], v['体积'], v['中心'][2])
              for v in 几何['体']}
    状态.多孔体 = st.multiselect('多孔体（管束），没有就留空', list(体标签), default=状态.多孔体,
                                format_func=体标签.get, key='多孔体_%d' % 状态.边界版本)
    标注文件 = st.file_uploader('或者上传 标注.json 覆盖', type=['json'], key='上传标注')
    if 标注文件 is not None and st.button('套用 标注.json', key='套用标注'):
        try:
            标注 = json.loads(标注文件.getvalue().decode('utf-8'))
            for t, 类 in 标注.get('面', {}).items():
                if int(t) not in 状态.面类 or 状态.自动面类[int(t)] == 'interface':
                    raise ValueError('面 %s 不存在或是内部面' % t)
                if 类 not in ui_tools.可选面类型:
                    raise ValueError('面 %s 的类型 %r 不支持' % (t, 类))
                状态.面类[int(t)] = 类
            if '多孔体' in 标注:
                状态.多孔体 = [int(v) for v in 标注['多孔体']]
            状态.边界版本 += 1
            st.rerun()
        except Exception as 错:
            st.error('标注.json 用不了: %s' % 错)

with 左:
    try:
        if '三角网格' not in 状态:
            with gmsh锁():
                状态.三角网格 = ui_tools.面三角网格(状态.step路径)
        st.plotly_chart(ui_tools.预览图(状态.三角网格, 几何, 状态.面类, 状态.多孔体), width='stretch')
    except ImportError:
        st.info('装上 plotly（pip install plotly）后这里会显示三维预览。')
    except Exception as 错:
        st.warning('三维预览画不出来: %s' % 错)

计数 = {}
for f in 外表面:
    计数[状态.面类[f['编号']]] = 计数.get(状态.面类[f['编号']], 0) + 1
st.markdown('**当前面类型**: %s；多孔体: %s'
            % ('，'.join('%s %d 个' % (k, 计数[k]) for k in sorted(计数)), 状态.多孔体 or '无'))
边界问题 = [名 for 名 in ('inlet', 'outlet') if not 计数.get(名)]
if 边界问题:
    st.error('还没有 %s 面，第 3 步会不通过。请在表里指定。' % ' 和 '.join(边界问题))

# ================================================================ 3 工况
st.header('3  工况')
工况 = 状态.工况
v = 状态.工况版本
上传列, 下载列 = st.columns(2)
工况文件 = 上传列.file_uploader('上传 工况.json（缺的项用当前值，多余的项忽略）', type=['json'], key='上传工况')
if 工况文件 is not None and 上传列.button('读入 工况.json', key='读入工况'):
    try:
        新, 不认识 = ui_tools.合并工况(工况, json.loads(工况文件.getvalue().decode('utf-8')))
        状态.工况, 状态.工况版本 = 新, v + 1
        if 不认识:
            st.warning('忽略了不认识的项: %s' % '、'.join(不认识))
        st.rerun()
    except Exception as 错:
        st.error('工况.json 用不了: %s' % 错)


def 数(组, 键, 标签, 列, 最小=None, 格式='%.6g', 帮助=None):
    组[键] = 列.number_input(标签, value=float(组[键]), min_value=最小, format=格式, help=帮助,
                            key='工况_%d_%s' % (v, 标签))


c = st.columns(3)
数(工况['入口'], '速度', '入口速度 (m/s，沿入口面法向)', c[0], 0.0)
数(工况['入口'], '温度', '入口温度 (K)', c[1], 0.0)
数(工况['出口'], '表压', '出口表压 (Pa)', c[2])
with st.expander('物性和壁面', expanded=False):
    c = st.columns(3)
    数(工况['物性'], '工作压力', '工作压力 (Pa)', c[0], 0.0)
    数(工况['物性'], '摩尔质量', '摩尔质量 (kg/mol)', c[1], 0.0)
    数(工况['物性'], '比热', '比热 (J/(kg·K))', c[2], 0.0)
    c = st.columns(3)
    数(工况['物性'], '湍流粘度', '湍流粘度 (Pa·s，常数)', c[0], 0.0,
       帮助='层流求解加一个常数湍流粘度；它决定有效粘度和导热')
    数(工况['物性'], '湍流普朗特数', '湍流普朗特数', c[1], 0.0)
    工况['壁面'] = c[2].selectbox('壁面', ('滑移', '无滑移'), index=('滑移', '无滑移').index(工况['壁面']),
                                key='工况_%d_壁面' % v)
if 状态.多孔体:
    st.subheader('多孔区（管束）')
    孔 = 工况.setdefault('多孔', ui_tools.默认工况()['多孔'])
    c = st.columns(3)
    数(孔, '粘性阻力系数', '粘性阻力系数 (1/m²)', c[0], 0.0)
    数(孔, '惯性阻力系数', '惯性阻力系数 (1/m)', c[1], 0.0)
    数(孔, '横向倍数', '横向倍数（横向阻力是流向的几倍）', c[2], 0.0)
    换 = 孔['换热']
    换['方式'] = st.radio('换热方式', ('温差', '定热源'), index=('温差', '定热源').index(换['方式']),
                        horizontal=True, key='工况_%d_换热方式' % v)
    c = st.columns(3)
    if 换['方式'] == '温差':
        换.setdefault('h_av', 500.0), 换.setdefault('壁温顶部', 500.0), 换.setdefault('壁温底部', 450.0)
        数(换, 'h_av', '体积换热系数 h_av (W/(m³·K))', c[0], 0.0)
        数(换, '壁温顶部', '壁温: 多孔区流入端 (K)', c[1], 0.0)
        数(换, '壁温底部', '壁温: 多孔区流出端 (K)', c[2], 0.0)
    else:
        换.setdefault('体积热源', -1e5)
        数(换, '体积热源', '体积热源 (W/m³，吸热为负)', c[0])
工况最终 = ui_tools.整理工况(工况, bool(状态.多孔体))
下载列.download_button('下载当前 工况.json', ui_tools.转json(工况最终), file_name='工况.json',
                       mime='application/json', key='下载工况')

# ================================================================ 4 计算设置
st.header('4  计算设置')
c = st.columns(4)
默认尺寸 = 几何['对角线'] / 30
网格尺寸 = c[0].number_input('网格尺寸 (m)，0 表示自动', value=0.0, min_value=0.0, format='%.4g',
                            help='自动取包围盒对角线的 1/30，当前为 %.4g m' % 默认尺寸, key='网格尺寸')
步数 = c[1].number_input('求解步数', value=默认步数, min_value=50, step=100, key='步数')
单元上限 = c[2].number_input('单元数上限', value=min(auto_run.默认单元上限, 限额['单元']), min_value=10000,
                            max_value=限额['单元'], step=50000,
                            help='自动加密或网格无关性检查估计超过这个数就停下', key='单元上限')
最多返工 = c[3].number_input('每层网格最多自动返工次数', value=auto_run.默认最多返工, min_value=1, max_value=10,
                            key='最多返工')
c = st.columns(4)
网格检查 = c[0].checkbox('网格无关性检查', value=True, key='网格检查',
                        help='通过后把网格加密一层再算，压差变化不超过 3%、吸热量不超过 2% 才算结果与网格无关')
对照一维 = c[1].checkbox('对照一维估算', value=False, key='对照一维',
                        help='只适用于等截面、来流均匀的多孔区（如 S1）。变截面、转弯入口不要勾')
if 固定运行方式:
    运行方式 = 固定运行方式
    c[2].markdown('运行方式：%s（服务器设定）' % 运行方式)
else:
    运行方式 = c[2].radio('运行方式', ('容器', '本机'), horizontal=True, key='运行方式',
                          help='容器: 用 docker 运行 %s；本机: 直接用本机装的 OpenFOAM' % 默认镜像)
时限 = c[3].number_input('单次求解时限 (秒)', value=min(默认时限, 限额['时限']), min_value=60, max_value=限额['时限'],
                        step=600, key='时限', help='超时会自动加倍重试，最多到 7200 秒')

if not 固定运行方式 and st.button('检查运行环境', key='检查环境'):
    if 运行方式 == '容器':
        if shutil.which('docker') is None:
            st.error('找不到 docker。请先安装 Docker Desktop 并启动，或改用“本机”。')
        else:
            信息 = subprocess.run(['docker', 'info'], capture_output=True, text=True)
            if 信息.returncode != 0:
                st.error('docker 没有在运行: %s' % (信息.stderr.strip()[:300]))
            else:
                镜像 = subprocess.run(['docker', 'image', 'inspect', 默认镜像], capture_output=True, text=True)
                st.success('docker 可用。' + ('镜像 %s 已在本机。' % 默认镜像 if 镜像.returncode == 0 else
                                            '镜像 %s 还没下载，第一次计算会先自动下载（约 1 GB）。' % 默认镜像))
    else:
        有 = [p for p in ('/usr/lib/openfoam', '/opt') if os.path.isdir(p) and
             any(n.startswith('openfoam') for n in os.listdir(p))]
        st.success('找到 OpenFOAM: %s' % 有) if 有 else st.error('没有在 /usr/lib/openfoam 或 /opt 下找到 OpenFOAM')

# ================================================================ 5 提交
st.header('5  提交计算')
标注 = {'面': {str(f['编号']): 状态.面类[f['编号']] for f in 外表面}, '多孔体': sorted(状态.多孔体)}
设置 = {'网格尺寸': 网格尺寸 or None, '步数': int(步数), '对照一维': 对照一维, '网格检查': 网格检查,
        '运行方式': 运行方式, '时限': int(时限), '单元上限': int(单元上限), '最多返工': int(最多返工)}
说明 = st.text_input('任务说明（可选，方便以后找）', value=状态.get('几何说明', ''), key='任务说明')
st.caption('提交后在后台排队计算，网页可以关掉；到“我的任务”看进度和结果。')
if st.button('提交计算', type='primary', disabled=bool(边界问题), key='提交计算'):
    当前 = 指纹(状态.step路径, 标注, 工况最终, 设置)
    if 状态.get('上次提交') == 当前 and st.session_state.get('上次任务'):
        st.warning('同样的输入刚刚已经提交过（任务 %s），没有重复提交。' % 状态.上次任务)
    else:
        try:
            号 = job_queue.提交(状态.用户, 状态.step路径, 工况最终, 标注, 设置, 说明=说明, 每人上限=限额['每人任务'])
        except job_queue.任务错误 as 错:
            st.error(str(错))
        else:
            状态.上次提交, 状态.上次任务 = 当前, 号
            状态.选中任务, 状态.跳转任务 = 号, True
            st.rerun()
