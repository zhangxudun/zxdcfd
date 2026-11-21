# -*- coding: utf-8 -*-
'''端到端测试: 第 1 步生成的 S1–S3 一路交给第 2 至 8 步，中间不经过人工。

链路: generate_step -> generate_boundary -> validate_boundary -> generate_mesh
      -> generate_case -> run_case -> check_result -> judge_result。
测试点 E2E-01 至 E2E-06 对应 test_E2E_01 至 test_E2E_06:
    E2E-01 至 E2E-03  第 1–2 步，S1–S3 各跑一次，不需要 OpenFOAM。
    E2E-04            第 3–5 步，S1–S3 各跑一次，不需要 OpenFOAM。
    E2E-05            第 6–8 步，S1–S3 不对照一维，四项检查都通过，结论是“通过”。
    E2E-06            S1–S3 再做一次一维对照，六项都通过，结论仍是“通过”。压差取多孔区前后截面，变截面、转弯都能比；
                      S3 在默认网格上出口烟温偏差超过 3%（一阶离散误差），用 一维对照网格比例 加密一层再算。
第 2 步的标注由 generate_structure.结构标注 按结构卡片给出: S1、S2 没有（返回 None，全自动识别），
S3 给出入口面和顶面，其余仍自动识别。“中间不经过人工”指不手写面编号。
E2E-05、E2E-06 的依据是 docs/07-检查结果.md 和 docs/08-判定与回退.md 的第 4 节。
它们用容器真的运行 OpenFOAM，步数用第 5 步的默认值，比较慢；每种结构只求解一次，各测试共用。
不能运行 docker 时这两个测试跳过（跳过不算通过），CI 里不允许跳过。
E2E-04 至 E2E-06 每算完一步就把实际结果记进端到端报告（tests/e2e_report.py），测完写到 E2E_REPORT 指定的文件，
没有指定时写到仓库根目录的 e2e-report.txt。报告里的数值就是这次测试判定用的那一次，不另外再算。
预期值按 docs/01-01-生成三维结构.md 结构卡片里的公式由参数算出。面和体按位置去找，不写死编号。
运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
'''
import collections
import json
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_structure import generate_step, 默认参数, 结构清单, 结构标注  # noqa: E402
from generate_boundary import generate_boundary, 标注文件名  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402
from generate_mesh import generate_mesh  # noqa: E402
from generate_case import generate_case, 默认步数  # noqa: E402
from run_case import run_case  # noqa: E402
from check_result import check_result  # noqa: E402
from judge_result import judge_result  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2e_report as 报  # noqa: E402

# E2E-06 里这些结构先加密一层网格再做一维对照。S3 默认网格（0.765 m）上出口烟温比一维估算高 3.29 K，
# 超过估算温降的 3%（3.04 K）；网格尺寸 0.765、0.547、0.39 m 时偏差为 3.29、约 2.5、1.86 K，随加密单调减小，
# 外推到零网格尺寸约 0.15 K，所以是一阶迎风的离散误差，不是一维估算不适用。取 1/1.4（与自动运行的加密比相同）。
一维对照网格比例 = {'S3': 1 / 1.4}
相对容差 = 1e-6
绝对容差 = 1e-3
网格容差 = 0.01          # 第 4 步的通过条件: 网格体积与几何体积的相对偏差不超过 1%
全部结构 = list(结构清单)
基本检查 = ['流量守恒', '入口流量', '监测值稳定', '残差收敛']
一维检查 = ['压差与一维解一致', '出口烟温与一维解一致']
边界类型 = {'inlet': 'patch', 'outlet': 'patch', 'symmetry': 'symmetry', 'walls': 'wall'}
算例文件 = ('Allrun', 'mesh.msh', '算例信息.json', '0/U', '0/p', '0/T',
            'system/controlDict', 'system/fvOptions', 'constant/transportProperties')


# ---------------------------------------------------------------- 工况与预期

def 多孔工况():
    return {'入口': {'速度': 5.0, '温度': 600.0},
            '出口': {'表压': 0.0},
            '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                     '湍流粘度': 1e-3, '湍流普朗特数': 0.85},
            '壁面': '滑移',
            '多孔': {'粘性阻力系数': 1e4, '惯性阻力系数': 10.0, '横向倍数': 100.0,
                     '换热': {'方式': '温差', 'h_av': 500.0,
                              '壁温顶部': 500.0, '壁温底部': 450.0}}}


def 预期_S1(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    return {'计数': {'inlet': 1, 'outlet': 1, 'interface': 2, 'walls': 12},
            '入口': ((0, b / 2, h + li), 2 * a * b),
            '出口': ((0, b / 2, -lo), 2 * a * b),
            '管束': ((0, b / 2, h / 2), 2 * a * b * h),
            '共用面': ((0.0, h), 2 * a * b),
            '对称面': (0, 0.0)}


def 预期_S2(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    ai, bi = p['入口半宽x'], p['入口深度y']
    ao, bo = p['出口半宽x'], p['出口深度y']
    return {'计数': {'inlet': 1, 'outlet': 1, 'interface': 2, 'symmetry': 3, 'walls': 9},
            '入口': ((ai / 2, b / 2, h + li), ai * bi),
            '出口': ((ao / 2, b / 2, -lo), ao * bo),
            '管束': ((a / 2, b / 2, h / 2), a * b * h),
            '共用面': ((0.0, h), a * b),
            '对称面': (3, (b + bi) / 2 * li + b * h + (b + bo) / 2 * lo)}


def 预期_S3(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    lh, hd = p['入口烟道长度'], p['入口烟道高度']
    return {'计数': {'inlet': 1, 'outlet': 1, 'interface': 2, 'walls': 14},
            '入口': ((a + lh, b / 2, h + li - hd / 2), b * hd),
            '出口': ((0, b / 2, -lo), 2 * a * b),
            '管束': ((0, b / 2, h / 2), 2 * a * b * h),
            '共用面': ((0.0, h), 2 * a * b),
            '对称面': (0, 0.0)}


预期 = {'S1': 预期_S1, 'S2': 预期_S2, 'S3': 预期_S3}


# ---------------------------------------------------------------- 工具

def 近似(实际, 预期值, 容差=相对容差):
    return abs(实际 - 预期值) <= 容差 * abs(预期值)


def 位置一致(实际, 预期值):
    return all(abs(实际[i] - 预期值[i]) <= 绝对容差 for i in range(3))


def 按类型(结果, 类):
    return [f for f in 结果['几何']['面'] if 结果['面'][f['编号']] == 类]


def 读文件(路径):
    with open(路径, encoding='utf-8') as f:
        return json.load(f)


def 报告(结果):
    return json.dumps(结果, ensure_ascii=False, indent=1, default=str)


def 不记(行):
    '''不需要报告时的记录函数（E2E-01 至 E2E-03），什么都不记。'''


def 跑两步(结构号, 参数, 目录, 记=不记):
    '''第 1 步生成 STEP，第 2 步直接读这个 STEP（带上结构卡片给的标注）。返回 (STEP 路径, 第 2 步结果)。'''
    模型 = generate_step(结构号, 参数, 目录 / 'geo', 'model')
    assert os.path.isfile(模型) and os.path.getsize(模型) > 0
    标注 = 结构标注(结构号, 参数, 模型)
    报.记第1步(记, 模型, 标注)
    边界 = generate_boundary(模型, 多孔工况(), 目录 / 'bc', 标注)
    报.记第2步(记, 边界)
    return 模型, 边界


def 跑前五步(结构号, 目录, 记=不记, 网格比例=None):
    '''默认参数走第 1–5 步，每一步直接用上一步的返回值。步数用默认值；网格尺寸用默认值，
    给了 网格比例 时取默认尺寸乘这个比例。

    记 是端到端报告的记录函数，每算完一步记一次。'''
    参数 = 默认参数(结构号)
    报.记输入(记, 结构号, 参数, 多孔工况())
    模型, 边界 = 跑两步(结构号, 参数, 目录, 记)
    校验 = validate_boundary(边界)
    报.记第3步(记, 校验)
    网格 = generate_mesh(模型, 边界, 目录 / 'mesh')
    if 网格比例 is not None:
        网格 = generate_mesh(模型, 边界, 目录 / 'mesh', 网格['网格尺寸'] * 网格比例)
    报.记第4步(记, 网格)
    算例 = generate_case(模型, 边界, 网格, 目录 / '算例')
    报.记第5步(记, 算例)
    return {'参数': 参数, '模型': 模型, '边界': 边界, '校验': 校验, '网格': 网格, '算例': 算例}


def docker可用():
    if shutil.which('docker') is None:
        return False
    return subprocess.run(['docker', 'info'], stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0


def 只算一次(已有, 键, 算):
    '''失败也记下来，同一种结构不重复划网格、不重复求解。'''
    if 键 not in 已有:
        try:
            已有[键] = 算()
        except Exception as 错:
            已有[键] = 错
    if isinstance(已有[键], Exception):
        raise 已有[键]
    return 已有[键]


def 记下中断(记, 算):
    '''算() 报错时先把中断和第 8 步的返工结论记进报告，再照原样抛出。'''
    try:
        return 算()
    except Exception as 错:
        报.记中断(记, 错)
        raise


@pytest.fixture(scope='module')
def 端到端报告():
    '''整个模块共用一份报告，模块里的测试都跑完后写出。'''
    报告 = 报.端到端报告()
    yield 报告
    报告.写出()


@pytest.fixture(scope='module')
def 前五步(tmp_path_factory, 端到端报告):
    '''前五步(结构号) 返回 跑前五步 的结果。'''
    已有 = {}

    def 取(结构号):
        记 = 端到端报告.记录(结构号)
        return 只算一次(已有, 结构号, lambda: 记下中断(
            记, lambda: dict(跑前五步(结构号, tmp_path_factory.mktemp(结构号), 记), 记=记)))
    return 取


@pytest.fixture(scope='module')
def 加密前五步(tmp_path_factory, 端到端报告):
    '''加密前五步(结构号) 按 一维对照网格比例 加密网格后跑前五步，报告里单独成段。'''
    已有 = {}

    def 取(结构号):
        段 = '%s（一维对照，网格 × %.4g）' % (结构号, 一维对照网格比例[结构号])
        记 = 端到端报告.记录(段)
        return 只算一次(已有, 结构号, lambda: 记下中断(记, lambda: dict(
            跑前五步(结构号, tmp_path_factory.mktemp(结构号 + '-fine'), 记, 一维对照网格比例[结构号]), 记=记)))
    return 取


@pytest.fixture(scope='module')
def 全流程(前五步):
    '''全流程(结构号) 在 前五步 的结果上加一项 求解（第 6 步的返回值）。'''
    已有, 记过跳过 = {}, set()

    def 求解(数据):
        结果 = run_case(数据['算例'])
        报.记第6步(数据['记'], 结果)
        return dict(数据, 求解=结果)

    def 取(结构号):
        数据 = 前五步(结构号)
        if not docker可用():
            if 结构号 not in 记过跳过:
                记过跳过.add(结构号)
                数据['记']('[第 6 步] 跳过: 不能运行 docker')
            if os.environ.get('CI'):
                pytest.fail('CI 里不允许跳过 OF 类测试: 不能运行 docker')
            pytest.skip('不能运行 docker，OF 类测试跳过（跳过不算通过）')
        return 只算一次(已有, 结构号, lambda: 记下中断(数据['记'], lambda: 求解(数据)))
    return 取


def 检查识别(结果, 卡):
    几何 = 结果['几何']
    assert dict(collections.Counter(结果['面'].values())) == 卡['计数']

    for 类, 名 in (('inlet', '入口'), ('outlet', '出口')):
        面 = 按类型(结果, 类)
        assert len(面) == 1, 名
        中心, 面积 = 卡[名]
        assert 位置一致(面[0]['中心'], 中心), 名
        assert 近似(面[0]['面积'], 面积), 名
        assert len(面[0]['相邻体']) == 1, 名

    assert len(几何['体']) == 3
    assert len(结果['多孔体']) == 1
    管束 = [v for v in 几何['体'] if v['编号'] == 结果['多孔体'][0]]
    assert len(管束) == 1
    中心, 体积 = 卡['管束']
    assert 位置一致(管束[0]['中心'], 中心)
    assert 近似(管束[0]['体积'], 体积)

    共用z, 共用面积 = 卡['共用面']
    内部 = sorted(按类型(结果, 'interface'), key=lambda f: f['中心'][2])
    assert len(内部) == len(共用z)
    for 面, z in zip(内部, sorted(共用z)):
        assert abs(面['中心'][2] - z) <= 绝对容差
        assert 近似(面['面积'], 共用面积)
        assert len(面['相邻体']) == 2 and 结果['多孔体'][0] in 面['相邻体']

    个数, 面积合计 = 卡['对称面']
    对称 = 按类型(结果, 'symmetry')
    assert len(对称) == 个数
    if 个数:
        assert all(abs(f['中心'][0]) <= 绝对容差 for f in 对称)
        assert 近似(sum(f['面积'] for f in 对称), 面积合计)


def 检查结论(数据, 对照一维):
    '''第 7 步全部通过，第 8 步结论是“通过”。返回第 7 步的结果。'''
    求解, 算例, 记 = 数据['求解'], 数据['算例'], 数据['记']
    结果 = 记下中断(记, lambda: check_result(求解, 对照一维=对照一维))
    报.记第7步(记, 结果, 对照一维)
    结论 = 记下中断(记, lambda: judge_result(结果))
    报.记结论(记, 结论)
    assert [项['项目'] for 项 in 结果['检查']] == 基本检查 + (一维检查 if 对照一维 else []), 报告(结果)
    assert all(项['通过'] is True for 项 in 结果['检查']), 报告(结果)
    assert 结果['通过'] is True, 报告(结果)
    assert 结果['迭代步数'] == 求解['完成步数']
    assert 结果['多孔压降估算'] == 算例['多孔压降估算']
    assert 结果['出口温度估算'] == 算例['出口温度估算']
    assert 读文件(结果['结果路径'])['通过'] is True

    assert 结论 == {'结论': '通过', '回到步骤': [], '对照一维': 对照一维, '原因': []}, 报告(结论)
    assert judge_result(求解['算例目录']) == 结论
    return 结果


# ---------------------------------------------------------------- 测试

@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_01(结构号, tmp_path):
    '''默认参数: 只带结构卡片给的标注（S1、S2 没有），面类型、入口、出口、多孔体、对称面全部认对，工况通过。'''
    参数 = 默认参数(结构号)
    _, 结果 = 跑两步(结构号, 参数, tmp_path)
    检查识别(结果, 预期[结构号](参数))
    assert 结果['工况'] == 多孔工况()


@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_02(结构号, tmp_path):
    '''标注.json 只含外表面且都有类型；拿它当标注再跑一次，结果不变。'''
    模型, 结果 = 跑两步(结构号, 默认参数(结构号), tmp_path)
    assert os.path.basename(结果['标注路径']) == 标注文件名
    文件 = 读文件(结果['标注路径'])
    外表面 = {t: 类 for t, 类 in 结果['面'].items() if 类 != 'interface'}
    卡 = 预期[结构号](默认参数(结构号))
    assert len(文件['面']) == sum(卡['计数'].values()) - 卡['计数']['interface']
    assert {int(k): v for k, v in 文件['面'].items()} == 外表面
    assert 文件['多孔体'] == 结果['多孔体']

    再跑 = generate_boundary(模型, 多孔工况(), tmp_path / 'bc2', 结果['标注路径'])
    assert 再跑['面'] == 结果['面']
    assert 再跑['多孔体'] == 结果['多孔体']
    assert 读文件(再跑['标注路径']) == 文件


@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_03(结构号, tmp_path):
    '''改参数（管束高度 3.2 改为 5.0）后，识别结果跟着参数变。'''
    参数 = dict(默认参数(结构号), 管束高度=5.0)
    _, 结果 = 跑两步(结构号, 参数, tmp_path)
    卡 = 预期[结构号](参数)
    检查识别(结果, 卡)
    assert abs(按类型(结果, 'inlet')[0]['中心'][2] - 卡['入口'][0][2]) <= 绝对容差
    assert abs(卡['入口'][0][2] - {'S1': 8.0, 'S2': 8.0, 'S3': 9.0}[结构号]) <= 绝对容差
    assert 近似(卡['管束'][1], {'S1': 480.0, 'S2': 240.0, 'S3': 480.0}[结构号])


@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_04(结构号, 前五步):
    '''第 2 步的结果直接交给第 3、4、5 步: 校验通过，网格和算例里的分组、面积、体积、流量都对得上。'''
    数据 = 前五步(结构号)
    参数, 边界, 校验, 网格, 算例 = (数据[键] for 键 in ('参数', '边界', '校验', '网格', '算例'))
    卡 = 预期[结构号](参数)
    工况 = 多孔工况()
    外表面类型 = {类 for 类 in 卡['计数'] if 类 != 'interface'}

    # 第 3 步: 入口、出口就是第 2 步认出来的那两个面
    assert 校验 == {类: [按类型(边界, 类)[0]['编号']] for 类 in ('inlet', 'outlet')}

    # 第 4 步: 分组名与第 2 步的面类型、多孔体一致，体积和面积与结构卡片一致
    assert os.path.isfile(网格['网格路径']) and os.path.getsize(网格['网格路径']) > 0
    多孔体 = 边界['多孔体']
    其余体 = sorted(v['编号'] for v in 边界['几何']['体'] if v['编号'] not in 多孔体)
    assert 网格['多孔体'] == 多孔体
    assert 网格['体分组'] == {'porous_zone': 多孔体, 'fluid_zone': 其余体}
    assert set(网格['边界']) == 外表面类型
    assert 网格['单元数'] > 0 and 网格['相对偏差'] <= 网格容差
    assert 近似(网格['体'][多孔体[0]]['网格体积'], 卡['管束'][1], 网格容差)
    assert 近似(网格['边界']['inlet']['面积'], 卡['入口'][1], 网格容差)
    assert 近似(网格['边界']['outlet']['面积'], 卡['出口'][1], 网格容差)

    # 第 5 步: 算例信息里的数与结构参数、工况一致，一维估算值算得出来
    assert 算例['边界'] == {类: 边界类型[类] for 类 in 外表面类型}
    assert 算例['多孔区'] == 'porous_zone' and 算例['步数'] == 默认步数
    assert 近似(算例['入口面积'], 卡['入口'][1])
    流量 = 工况['入口']['速度'] * 卡['入口'][1]
    assert 近似(算例['理论入口体积流量'], 流量)
    assert 近似(算例['多孔区长度'], 参数['管束高度'])
    assert 近似(算例['多孔区流速'], 流量 / 卡['共用面'][1])
    起, 终 = 算例['多孔范围']                    # 流向是 -z，沿流向的坐标是 -z
    assert abs(起 + 参数['管束高度']) <= 绝对容差 and abs(终) <= 绝对容差
    assert 算例['多孔压降估算'] > 0
    assert 工况['多孔']['换热']['壁温底部'] < 算例['出口温度估算'] < 工况['入口']['温度']
    for 名 in 算例文件:
        assert os.path.isfile(os.path.join(算例['算例目录'], *名.split('/'))), 名
    assert 读文件(算例['算例信息路径'])['理论入口体积流量'] == 算例['理论入口体积流量']


@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_05(结构号, 全流程):
    '''第 6–8 步: 求解跑完，不对照一维时四项检查都通过，结论是“通过”。'''
    数据 = 全流程(结构号)
    求解, 算例 = 数据['求解'], 数据['算例']
    assert 求解['边界'] == 算例['边界']
    assert 求解['步数'] == 默认步数
    assert 求解['完成步数'] == 默认步数 or (求解['提前收敛'] and 求解['完成步数'] < 默认步数)
    assert sorted(求解['监测数据']) == sorted(算例['监测量'])

    结果 = 检查结论(数据, 对照一维=False)
    工况 = 多孔工况()
    assert 工况['多孔']['换热']['壁温底部'] < 结果['出口温度'] < 工况['入口']['温度'], 报告(结果)
    assert 结果['吸热量'] > 0, 报告(结果)


@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_06(结构号, 全流程, 加密前五步):
    '''再做一次一维对照: 六项都通过，结论仍是“通过”。压差取多孔区前后截面（porousInP − porousOutP），
    不含进出口动压变化和转弯损失。一维对照网格比例 里的结构先加密网格再算。'''
    if 结构号 in 一维对照网格比例:
        数据 = 加密前五步(结构号)
        求解 = 记下中断(数据['记'], lambda: run_case(数据['算例']))
        报.记第6步(数据['记'], 求解)
        数据 = dict(数据, 求解=求解)
    else:
        数据 = 全流程(结构号)
    结果 = 检查结论(数据, 对照一维=True)
    assert 结果['多孔区压差'] is not None, 报告(结果)
