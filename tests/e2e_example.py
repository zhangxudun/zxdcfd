# -*- coding: utf-8 -*-
'''端到端示例: 用默认参数把清单里的结构各跑一遍第 1–8 步，打印每一步的实际结果。

这不是测试，不做断言，只是把 tests/test_end_to_end.py 走的那条链路的中间结果摆出来。
工况用 test_end_to_end.py 里的 多孔工况()；网格尺寸、步数都用默认值。
S1、S2 不加标注。S3 的入口在侧面，第 2 步自动识别会认错，标注由 generate_structure.结构标注 按结构卡片给出。
第 6 步用容器真的运行 OpenFOAM，比较慢。第 7 步 S1 对照一维；S2 是变截面、S3 来流不均匀，不对照。
第 6 步之后另外从 log.simpleFoam 里读出各方程的残差打印出来，方便看收敛过程；第 7 步的“残差收敛”用的是同一个窗口。
第 3–7 步有哪一步报错（比如本机找不到 docker），就把异常交给第 8 步，打印返工结论，接着跑下一种结构。
结论是返工时退出码仍是 0；通不通过由 test_end_to_end.py 判断。
CI 会把输出发布到 ci-report 分支的 e2e-example.txt（PR 发布到 ci-report-pr 分支）。

用法（在仓库根目录）:
    python tests/e2e_example.py [输出目录]
输出目录默认是 e2e-example，每种结构一个子目录，里面是 model.step、标注.json、mesh.msh 和 算例/。
'''
import collections
import json
import math
import os
import re
import sys

_这里 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _这里)
sys.path.insert(0, os.path.dirname(_这里))
from generate_structure import generate_step, 默认参数, 结构清单, 结构标注  # noqa: E402
from generate_boundary import generate_boundary  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402
from generate_mesh import generate_mesh  # noqa: E402
from generate_case import generate_case  # noqa: E402
from run_case import run_case  # noqa: E402
from check_result import check_result  # noqa: E402
from judge_result import judge_result, 判定错误  # noqa: E402
from test_end_to_end import 多孔工况  # noqa: E402

对照一维的结构 = ('S1',)      # S2 是变截面，一维对照是否适用未定，见 docs/07-检查结果.md 第 6 节
残差窗口 = 51                 # 与第 7 步判稳用的窗口一样长
_残差行 = re.compile(r'Solving for (\w+), Initial residual = ([-+.\deE]+)')


def _点(坐标):
    return '(%s)' % ', '.join('%g' % (round(v, 6) + 0.0) for v in 坐标)


def _值(v):
    return '无' if v is None else '%.6g' % v


def _文(v):
    return json.dumps(v, ensure_ascii=False)


def 读残差(日志路径):
    '''从 log.simpleFoam 读每个时间步里各方程第一次求解的初始残差，返回 {方程: [每步的值]}。

    一个时间步里 p 会解几次（非正交修正），只取第一次，与 residualControl 看的是同一个数。'''
    残差, 本步 = collections.OrderedDict(), set()
    with open(日志路径, encoding='utf-8', errors='replace') as f:
        for 行 in f:
            if 行.startswith('Time = '):
                本步 = set()
                continue
            找到 = _残差行.search(行)
            if not 找到 or 找到.group(1) in 本步:
                continue
            try:
                值 = float(找到.group(2))
            except ValueError:
                continue
            本步.add(找到.group(1))
            残差.setdefault(找到.group(1), []).append(值)
    return 残差


def 打印残差(日志路径):
    try:
        残差 = 读残差(日志路径)
    except OSError as 错:
        print('[第 6 步] 残差读不了: %s' % 错)
        return
    if not 残差:
        print('[第 6 步] %s 里没有读到残差' % 日志路径)
        return
    print('[第 6 步] 残差（每步第一次求解的初始残差；fvSolution 里提前结束的条件是 p、U 都低于 1e-06）:')
    for 名, 值列 in 残差.items():
        首, 末, 末段 = 值列[0], 值列[-1], 值列[-残差窗口:]
        下降 = '%.1f' % math.log10(首 / 末) if 首 > 0 and 末 > 0 else '算不出'
        print('  %-3s 共 %d 步  第 1 步 %.3g  最后一步 %.3g  最后 %d 步里最大 %.3g 最小 %.3g  全程最大 %.3g  下降 %s 个数量级'
              % (名, len(值列), 首, 末, len(末段), max(末段), min(末段), max(值列), 下降))


def 跑后五步(模型, 边界, 目录, 对照一维):
    '''第 3–7 步，每一步直接用上一步的返回值。返回第 7 步的结果。'''
    校验 = validate_boundary(边界)
    print('[第 3 步] inlet 面 %s，outlet 面 %s' % (校验['inlet'], 校验['outlet']))

    网格 = generate_mesh(模型, 边界, 目录)
    print('[第 4 步] 生成 %s: 四面体 %d 个，节点 %d 个，网格尺寸 %s m'
          % (网格['网格路径'], 网格['单元数'], 网格['节点数'], _值(网格['网格尺寸'])))
    print('[第 4 步] 网格体积 %s m³，几何体积 %s m³，相对偏差 %.3g%%'
          % (_值(网格['网格体积']), _值(网格['几何体积']), 100 * 网格['相对偏差']))
    for 体号, 体 in sorted(网格['体'].items()):
        print('  体 %d  四面体 %d 个  网格体积 %s m³  几何体积 %s m³  相对偏差 %.3g%%'
              % (体号, 体['单元数'], _值(体['网格体积']), _值(体['几何体积']), 100 * 体['相对偏差']))
    for 名, 组 in 网格['边界'].items():
        print('  %-9s  面 %s  网格面积 %s m²' % (名, 组['面'], _值(组['面积'])))
    print('[第 4 步] 体分组: %s' % _文(网格['体分组']))

    算例 = generate_case(模型, 边界, 网格, os.path.join(目录, '算例'))
    print('[第 5 步] 算例目录 %s，步数 %s，多孔区 %s' % (算例['算例目录'], 算例['步数'], 算例['多孔区']))
    print('[第 5 步] 边界: %s' % _文(算例['边界']))
    print('[第 5 步] 入口面积 %s m²，理论入口体积流量 %s m³/s，密度 %s kg/m³'
          % (_值(算例['入口面积']), _值(算例['理论入口体积流量']), _值(算例['密度'])))
    print('[第 5 步] 多孔区流速 %s m/s，多孔区长度 %s m，多孔压降估算 %s Pa，出口温度估算 %s K'
          % (_值(算例['多孔区流速']), _值(算例['多孔区长度']),
             _值(算例['多孔压降估算']), _值(算例['出口温度估算'])))

    求解 = run_case(算例)
    print('[第 6 步] 完成 %s 步（设定 %s 步）%s，监测量 %s'
          % (求解['完成步数'], 求解['步数'], '，提前收敛' if 求解['提前收敛'] else '',
             ' / '.join(sorted(求解['监测数据']))))
    print('[第 6 步] 转出的网格边界: %s' % _文(求解['边界']))
    打印残差(求解['日志']['simpleFoam'])

    结果 = check_result(求解, 对照一维=对照一维)
    print('[第 7 步] %s，写出 %s' % ('对照一维' if 对照一维 else '不对照一维', 结果['结果路径']))
    for 项 in 结果['检查']:
        print('  %s  %s: %s' % ('通过' if 项['通过'] else '未通过', 项['项目'], 项['说明']))
    print('[第 7 步] 迭代 %d 步，进出口压差 %s Pa，入口流量 %s m³/s，出口流量 %s m³/s，'
          '出口温度 %s K，吸热量 %s W'
          % (结果['迭代步数'], _值(结果['进出口压差']), _值(结果['入口流量']), _值(结果['出口流量']),
             _值(结果['出口温度']), _值(结果['吸热量'])))
    return 结果


def 打印结论(结论):
    for 项 in 结论['原因']:
        print('  %s，回到第 %s 步: %s。%s'
              % (项['现象'], '、'.join(str(n) for n in 项['回到步骤']), 项['处理'], 项['说明']))
    if 结论['结论'] == '通过':
        print('[第 8 步] 结论: 通过（%s）' % ('含一维对照' if 结论['对照一维'] else '未对照一维'))
    else:
        print('[第 8 步] 结论: 返工，回到第 %s 步' % '、'.join(str(n) for n in 结论['回到步骤']))


def 跑一个(结构号, 输出目录):
    参数 = 默认参数(结构号)
    目录 = os.path.join(输出目录, 结构号)
    对照一维 = 结构号 in 对照一维的结构
    print('=' * 20, 结构号, 结构清单[结构号]['名称'], '=' * 20)
    print('[输入] 结构参数: %s' % _文(参数))
    print('[输入] 工况: %s' % _文(多孔工况()))

    模型 = generate_step(结构号, 参数, 目录, 'model')
    print('[第 1 步] 生成 %s，%d 字节' % (模型, os.path.getsize(模型)))

    标注 = 结构标注(结构号, 参数, 模型)
    print('[输入] 标注: %s' % ('无' if 标注 is None else _文(标注)))

    结果 = generate_boundary(模型, 多孔工况(), 目录, 标注)
    几何 = 结果['几何']
    print('[第 2 步] 体 %d 个:' % len(几何['体']))
    for 体 in sorted(几何['体'], key=lambda v: v['编号']):
        print('  体 %d  体积 %.6g m³  中心 %s%s'
              % (体['编号'], 体['体积'], _点(体['中心']),
                 '  多孔体' if 体['编号'] in 结果['多孔体'] else ''))
    print('[第 2 步] 面 %d 个:' % len(几何['面']))
    for 面 in sorted(几何['面'], key=lambda f: f['编号']):
        print('  面 %2d  %-9s  面积 %.6g m²  中心 %s  相邻体 %s'
              % (面['编号'], 结果['面'][面['编号']], 面['面积'],
                 _点(面['中心']), 面['相邻体']))
    计数 = collections.Counter(结果['面'].values())
    print('[第 2 步] 面类型计数: %s'
          % ', '.join('%s %d' % 项 for 项 in sorted(计数.items())))
    print('[第 2 步] 多孔体: %s' % 结果['多孔体'])
    print('[第 2 步] 校验后的工况: %s' % _文(结果['工况']))
    print('[第 2 步] %s 内容:' % 结果['标注路径'])
    with open(结果['标注路径'], encoding='utf-8') as f:
        print(f.read().rstrip())

    try:
        结论 = judge_result(跑后五步(模型, 结果, 目录, 对照一维))
    except (ValueError, RuntimeError) as 错:
        try:
            结论 = judge_result(错)          # 第 2–7 步的报错，由第 8 步指出回到哪一步
        except 判定错误:
            raise 错 from None              # 不是流程里的报错，是程序出错，照原样抛出
        print('[中断] %s: %s' % (type(错).__name__, 错))
    打印结论(结论)
    print()


if __name__ == '__main__':
    输出目录 = sys.argv[1] if len(sys.argv) > 1 else 'e2e-example'
    for 结构号 in 结构清单:
        跑一个(结构号, 输出目录)
