# -*- coding: utf-8 -*-
'''docs/03-校验.md 的测试，测试点 T3-01 至 T3-14 对应 test_T3_01 至 test_T3_14。

运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
测试几何 A、C、E、G 在测试里用 gmsh 现建并写成 STEP，先经过第 2 步，再交给第 3 步。
面的编号按中心坐标去找，不写死。
'''
import copy
import os
import sys

import gmsh
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_boundary import generate_boundary, 识别几何, 标注文件名  # noqa: E402
from generate_structure import _写出_单位m  # noqa: E402
from validate_boundary import validate_boundary, 校验错误  # noqa: E402

绝对容差 = 1e-3
_轴 = {'x': 0, 'y': 1, 'z': 2}


# ---------------------------------------------------------------- 测试几何

def _长方体(occ):
    occ.addBox(0, 0, 0, 1, 2, 3)


def _叠放(occ):
    for z in range(3):
        occ.addBox(0, 0, z, 1, 2, 1)


def _弯头(occ):
    '''L 形截面在 y=0 平面上，沿 y 拉伸 1 m。'''
    角 = [(1, 0), (4, 0), (4, 1), (2, 1), (2, 3), (1, 3)]
    点 = [occ.addPoint(x, 0, z) for x, z in 角]
    边 = [occ.addLine(点[i], 点[(i + 1) % len(点)]) for i in range(len(点))]
    截面 = occ.addPlaneSurface([occ.addCurveLoop(边)])
    occ.extrude([(2, 截面)], 0, 1, 0)


def _横放圆柱(occ):
    '''轴线沿 y，过 x=1、z=0.5，半径 0.5 m，长 2 m。没有水平面。'''
    occ.addCylinder(1, 0, 0.5, 0, 2, 0, 0.5)


测试几何 = {'A': _长方体, 'C': _叠放, 'E': _弯头, 'G': _横放圆柱}


def 写STEP(路径, 构造):
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add('geo')
        构造(gmsh.model.occ)
        gmsh.model.occ.synchronize()
        _写出_单位m(路径)
    finally:
        gmsh.finalize()
    return 路径


@pytest.fixture(scope='module')
def 几何文件(tmp_path_factory):
    目录 = str(tmp_path_factory.mktemp('geo'))
    return {名: 写STEP(os.path.join(目录, 名 + '.step'), 构造) for 名, 构造 in 测试几何.items()}


# ---------------------------------------------------------------- 工况

def 基本工况():
    return {'入口': {'速度': 5.0, '温度': 600.0},
            '出口': {'表压': 0.0},
            '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                     '湍流粘度': 1e-3, '湍流普朗特数': 0.85},
            '壁面': '滑移'}


def 多孔工况():
    return dict(基本工况(), 多孔={
        '粘性阻力系数': 1e4, '惯性阻力系数': 10.0, '横向倍数': 100.0,
        '换热': {'方式': '温差', 'h_av': 500.0, '壁温顶部': 500.0, '壁温底部': 450.0}})


# ---------------------------------------------------------------- 工具

def 找面(几何, **坐标):
    '''按面中心的坐标找面编号，例如 找面(几何, z=3)。必须恰好命中一个。'''
    命中 = [f['编号'] for f in 几何['面']
            if all(abs(f['中心'][_轴[k]] - v) <= 绝对容差 for k, v in 坐标.items())]
    assert len(命中) == 1, '中心在 %r 的面应恰好一个，实际 %d 个' % (坐标, len(命中))
    return 命中[0]


def 第2步(几何文件, 名, 目录, 改=None, 工况=None):
    '''跑第 2 步。改 是 {面中心坐标: 类型}，例如 {('z', 3): 'walls'}，按坐标找到面再写成标注。'''
    标注 = None
    if 改:
        几何 = 识别几何(几何文件[名])
        标注 = {'面': {找面(几何, **{轴: 值}): 类 for (轴, 值), 类 in 改.items()}}
    return generate_boundary(几何文件[名], 工况 or 基本工况(), 目录, 标注)


def 目录快照(目录):
    return {n: open(os.path.join(目录, n), 'rb').read() for n in sorted(os.listdir(目录))}


def 报错信息(输入):
    with pytest.raises(校验错误) as e:
        validate_boundary(输入)
    return str(e.value)


去掉入口 = {('z', 3): 'walls'}
去掉出口 = {('z', 0): 'walls'}


# ---------------------------------------------------------------- 3.1 通过

def test_T3_01(几何文件, tmp_path):
    结果 = 第2步(几何文件, 'A', tmp_path)
    几何 = 结果['几何']
    assert validate_boundary(结果) == {'inlet': [找面(几何, z=3)], 'outlet': [找面(几何, z=0)]}


def test_T3_02(几何文件, tmp_path):
    结果 = 第2步(几何文件, 'C', tmp_path, 工况=多孔工况())
    几何 = 结果['几何']
    内部面 = {t for t, 类 in 结果['面'].items() if 类 == 'interface'}
    assert len(内部面) == 2
    统计 = validate_boundary(结果)
    assert 统计 == {'inlet': [找面(几何, z=3)], 'outlet': [找面(几何, z=0)]}
    assert not 内部面 & set(统计['inlet'] + 统计['outlet'])


def test_T3_03(几何文件, tmp_path):
    '''已知限制: 弯头的出口被猜在 z=0，是错的，但第 3 步只看有没有。规则改了要同步改这一条。'''
    结果 = 第2步(几何文件, 'E', tmp_path)
    几何 = 结果['几何']
    assert validate_boundary(结果) == {'inlet': [找面(几何, z=3)], 'outlet': [找面(几何, z=0)]}


def test_T3_04(几何文件, tmp_path):
    结果 = 第2步(几何文件, 'E', tmp_path, {('z', 0): 'walls', ('x', 4): 'outlet'})
    几何 = 结果['几何']
    assert validate_boundary(结果) == {'inlet': [找面(几何, z=3)], 'outlet': [找面(几何, x=4)]}


def test_T3_05(几何文件, tmp_path):
    结果 = 第2步(几何文件, 'A', tmp_path, {('y', 0): 'inlet'})
    几何 = 结果['几何']
    统计 = validate_boundary(结果)
    assert 统计['inlet'] == sorted([找面(几何, z=3), 找面(几何, y=0)])
    assert 统计['outlet'] == [找面(几何, z=0)]


def test_T3_06(几何文件, tmp_path):
    '''T3-09 报错后的修正路径: 回到第 2 步，在标注里指定面编号。'''
    结果 = 第2步(几何文件, 'G', tmp_path, {('y', 0): 'inlet', ('y', 2): 'outlet'})
    几何 = 结果['几何']
    assert validate_boundary(结果) == {'inlet': [找面(几何, y=0)], 'outlet': [找面(几何, y=2)]}


# ---------------------------------------------------------------- 3.2 不通过

def test_T3_07(几何文件, tmp_path):
    信息 = 报错信息(第2步(几何文件, 'A', tmp_path, 去掉入口))
    assert '没有 inlet 面' in 信息 and 'outlet' not in 信息


def test_T3_08(几何文件, tmp_path):
    信息 = 报错信息(第2步(几何文件, 'A', tmp_path, 去掉出口))
    assert '没有 outlet 面' in 信息 and 'inlet' not in 信息


def test_T3_09(几何文件, tmp_path):
    结果 = 第2步(几何文件, 'G', tmp_path)
    assert set(结果['面'].values()) == {'walls'}
    信息 = 报错信息(结果)
    assert 'inlet' in 信息 and 'outlet' in 信息


def test_T3_10(几何文件, tmp_path):
    信息 = 报错信息(第2步(几何文件, 'E', tmp_path, 去掉出口))
    assert '没有 outlet 面' in 信息 and 'inlet' not in 信息


# ---------------------------------------------------------------- 3.3 输入与产物

def test_T3_11(几何文件, tmp_path):
    通过 = 第2步(几何文件, 'A', tmp_path / 'ok')
    assert validate_boundary(通过['标注路径']) == validate_boundary(通过)
    assert validate_boundary(tmp_path / 'ok' / 标注文件名) == validate_boundary(通过)
    不通过 = 第2步(几何文件, 'A', tmp_path / 'bad', 去掉入口)
    assert 报错信息(不通过['标注路径']) == 报错信息(不通过)


def test_T3_12(tmp_path):
    不存在 = str(tmp_path / 'missing.json')
    坏文件 = str(tmp_path / 'bad.json')
    with open(坏文件, 'w', encoding='utf-8') as f:
        f.write('{"面": {')
    没有面 = str(tmp_path / 'noface.json')
    with open(没有面, 'w', encoding='utf-8') as f:
        f.write('{"多孔体": []}')
    assert 'missing.json' in 报错信息(不存在)
    assert 'JSON' in 报错信息(坏文件)
    assert '面' in 报错信息(没有面)
    assert '面' in 报错信息({'多孔体': []})


def test_T3_13(几何文件, tmp_path):
    for 名, 改, 应通过 in (('ok', None, True), ('bad', 去掉入口, False)):
        目录 = str(tmp_path / 名)
        结果 = 第2步(几何文件, 'A', 目录, 改)
        原结果, 之前 = copy.deepcopy(结果), 目录快照(目录)
        assert list(之前) == [标注文件名]
        for 输入 in (结果, 结果['标注路径']):
            if 应通过:
                validate_boundary(输入)
            else:
                报错信息(输入)
        assert 目录快照(目录) == 之前
        assert 结果 == 原结果


# ---------------------------------------------------------------- 3.4 可重复性

def test_T3_14(几何文件, tmp_path):
    通过 = 第2步(几何文件, 'A', tmp_path / 'ok')
    assert validate_boundary(通过) == validate_boundary(通过)
    不通过 = 第2步(几何文件, 'G', tmp_path / 'bad')
    assert 报错信息(不通过) == 报错信息(不通过)
