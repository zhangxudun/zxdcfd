# -*- coding: utf-8 -*-
'''网页界面 app.py 的测试（streamlit.testing.v1.AppTest，不开浏览器）。auto_run 用 tests/stub_auto_run.py 代替。

UI-01  单机模式（没设 ZXDCFD_USERS）不用登录，提示先给几何。
UI-02  选内置 S1 生成几何: 自动识别出 1 个 inlet、1 个 outlet，多孔体 1 个。
UI-03  提交计算: 任务进队列，面类型、工况、设置原样写进任务；后台算完后“我的任务”页显示结论和结果。
UI-04  内置 S3 生成几何: 套用结构卡片的标注，入口改到侧面，仍是 1 个 inlet。
UI-05  多人模式: 要登录；访问码错不能进；登录后只看得到自己的任务，管理员看得到全部。
UI-06  服务器限额: 固定运行方式时不显示选择；单元上限不能超过 ZXDCFD_MAX_CELLS；每人任务数到上限时提示。
UI-07  通用流体: 套用“水”的物性、外壁定壁温、均匀壁温后提交，任务里的工况是常数物性，能通过第 2 步的校验。
需要 pip install streamlit plotly gmsh。
'''
import os
import sys

import pytest

仓库 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, 仓库)
pytest.importorskip('streamlit')
from streamlit.testing.v1 import AppTest  # noqa: E402

import job_queue  # noqa: E402
import worker  # noqa: E402


@pytest.fixture
def 环境(tmp_path, monkeypatch):
    monkeypatch.setenv('ZXDCFD_WORKDIR', str(tmp_path / '运行目录'))
    monkeypatch.setenv('ZXDCFD_START_WORKER', '0')
    monkeypatch.setenv('ZXDCFD_AUTO_RUN', 'stub_auto_run:假运行')
    monkeypatch.setenv('PYTHONPATH', os.path.join(仓库, 'tests'))
    for 名 in ('ZXDCFD_USERS', 'ZXDCFD_USERS_FILE', 'ZXDCFD_ADMINS', 'ZXDCFD_RUN_MODE', 'ZXDCFD_MAX_CELLS',
               'ZXDCFD_MAX_JOBS_PER_USER', 'ZXDCFD_STUB_SLEEP', 'ZXDCFD_STUB_FAIL'):
        monkeypatch.delenv(名, raising=False)
    monkeypatch.setattr(worker, '轮询间隔', 0.2)
    return monkeypatch


def 新界面():
    return AppTest.from_file(os.path.join(仓库, 'app.py'), default_timeout=120)


def 文字(at, 类):
    return '\n'.join(str(e.value) for e in getattr(at, 类))


def 生成(at, 结构号):
    at.radio(key='几何来源').set_value('用内置结构生成').run()
    at.selectbox(key='结构号').set_value(结构号).run()
    at.button(key='生成几何').click().run()
    assert not at.exception, at.exception
    assert not at.error, 文字(at, 'error')


def 登录(at, 名, 码):
    at.run()
    at.text_input(key='登录名').input(名)
    at.text_input(key='登录码').input(码)
    at.button[0].click().run()
    return at


def test_UI_01_启动(环境):
    at = 新界面()
    at.run()
    assert not at.exception
    assert '先上传 STEP' in 文字(at, 'info')
    assert '单机模式' in 文字(at, 'markdown')


def test_UI_02_内置S1(环境):
    at = 新界面()
    at.run()
    生成(at, 'S1')
    assert '内置结构 S1' in 文字(at, 'success')
    面 = 文字(at, 'markdown')
    assert 'inlet 1 个' in 面 and 'outlet 1 个' in 面
    assert len(at.multiselect[0].value) == 1          # 多孔体


def test_UI_03_提交到结果(环境):
    at = 新界面()
    at.run()
    生成(at, 'S1')
    at.number_input(key='步数').set_value(1000)
    at.checkbox(key='网格检查').uncheck()
    at.button(key='提交计算').click().run()
    assert not at.exception, at.exception
    任务们 = job_queue.列出()
    assert len(任务们) == 1
    任务 = 任务们[0]
    assert 任务['用户'] == '本机' and 任务['状态'] == '排队'
    assert 任务['设置']['步数'] == 1000 and 任务['设置']['网格检查'] is False and 任务['设置']['网格尺寸'] is None
    assert at.radio(key='页面').value == '我的任务'
    assert '排队中' in 文字(at, 'markdown')

    worker.常驻(一次=True)
    assert job_queue.读(任务['号'])['状态'] == '完成'
    at.run()
    assert not at.exception, at.exception
    assert '替身说明' in 文字(at, 'success')
    assert any('22.9' in m.value for m in at.metric)
    import json
    with open(os.path.join(job_queue.任务目录(任务['号']), '输入', '标注.json'), encoding='utf-8') as f:
        标注 = json.load(f)
    assert sorted(标注['面'].values()).count('inlet') == 1 and len(标注['多孔体']) == 1


def test_UI_04_内置S3(环境):
    at = 新界面()
    at.run()
    生成(at, 'S3')
    面 = 文字(at, 'markdown')
    assert 'inlet 1 个' in 面 and 'outlet 1 个' in 面


def test_UI_05_多人登录与隔离(环境, tmp_path):
    环境.setenv('ZXDCFD_USERS', '甲:码甲,乙:码乙,管:码管')
    环境.setenv('ZXDCFD_ADMINS', '管')
    step = tmp_path / 'm.step'
    step.write_text('x')
    号甲 = job_queue.提交('甲', str(step), {}, {'面': {}, '多孔体': []}, {}, 说明='甲的任务')
    号乙 = job_queue.提交('乙', str(step), {}, {'面': {}, '多孔体': []}, {}, 说明='乙的任务')

    at = 登录(新界面(), '甲', '错码')
    assert '访问码不对' in 文字(at, 'error')
    assert not at.radio                               # 没进去

    at = 登录(新界面(), '甲', '码甲')
    assert not at.exception, at.exception
    at.radio(key='页面').set_value('我的任务').run()
    选项 = at.selectbox(key='选任务').options
    assert any(号甲 in o for o in 选项) and not any(号乙 in o for o in 选项)

    at = 登录(新界面(), '管', '码管')
    at.radio(key='页面').set_value('我的任务').run()
    选项 = at.selectbox(key='选任务').options
    assert any(号甲 in o for o in 选项) and any(号乙 in o for o in 选项)


def test_UI_06_服务器限额(环境):
    环境.setenv('ZXDCFD_RUN_MODE', '本机')
    环境.setenv('ZXDCFD_MAX_CELLS', '100000')
    环境.setenv('ZXDCFD_MAX_JOBS_PER_USER', '1')
    at = 新界面()
    at.run()
    生成(at, 'S1')
    assert '服务器设定' in 文字(at, 'markdown')
    assert at.number_input(key='单元上限').value == 100000
    at.button(key='提交计算').click().run()
    assert job_queue.列出()[0]['设置']['运行方式'] == '本机'
    at.radio(key='页面').set_value('新建计算').run()
    at.number_input(key='步数').set_value(1200)
    at.button(key='提交计算').click().run()
    assert '每人最多 1 个' in 文字(at, 'error')


def test_UI_07_通用流体(环境):
    from generate_boundary import 校验工况
    at = 新界面()
    at.run()
    生成(at, 'S1')
    at.selectbox(key='物性预设').set_value('水（20 °C 常数物性）').run()
    at.button(key='套用物性').click().run()
    assert not at.exception, at.exception
    v = at.session_state['工况版本']
    at.selectbox(key='工况_%d_壁面换热' % v).set_value('定壁温').run()
    at.radio(key='工况_%d_壁温写法' % v).set_value('均匀').run()
    at.checkbox(key='网格检查').uncheck()
    at.button(key='提交计算').click().run()
    assert not at.exception, at.exception
    工况 = job_queue.列出()[0]['工况'] if '工况' in job_queue.列出()[0] else None
    if 工况 is None:
        import json
        with open(os.path.join(job_queue.任务目录(job_queue.列出()[0]['号']), '输入', '工况.json'),
                  encoding='utf-8') as f:
            工况 = json.load(f)
    assert 工况['物性']['密度'] == 998.2 and '工作压力' not in 工况['物性']
    assert 工况['壁面换热']['方式'] == '定壁温' and '壁温' in 工况['多孔']['换热']
    校验工况(工况, True)
