# -*- coding: utf-8 -*-
'''docs/04-划分网格.md 的测试，测试点 T4-01 至 T4-19 对应 test_T4_01 至 test_T4_19。

运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
测试几何 A、C、E、G 在测试里用 gmsh 现建并写成 STEP，先经过第 2、3 步，再交给第 4 步。
面的编号按中心坐标去找，不写死。
'''
import copy
import math
import os
import sys

import gmsh
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_boundary import generate_boundary, 识别几何, 标注文件名, 几何错误, 标注错误  # noqa: E402
from generate_mesh import (generate_mesh, 网格错误, 网格文件名,  # noqa: E402
                           _四面体体积, _三角形面积)
from generate_structure import _写出_单位m  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402

绝对容差 = 1e-3
平面容差 = 1e-6     # 全是平面的几何: 体积、面积的相对误差
曲面容差 = 1e-2     # 有曲面的几何: 体积、面积的相对偏差
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
    '''轴线沿 y，过 x=1、z=0.5，半径 0.5 m，长 2 m。'''
    occ.addCylinder(1, 0, 0.5, 0, 2, 0, 0.5)


测试几何 = {'A': _长方体, 'C': _叠放, 'E': _弯头, 'G': _横放圆柱}
弯头标注 = {('z', 0): 'walls', ('x', 4): 'outlet'}
圆柱标注 = {('y', 0): 'inlet', ('y', 2): 'outlet'}


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
    命中 = [f['编号'] for f in 几何['面']
            if all(abs(f['中心'][_轴[k]] - v) <= 绝对容差 for k, v in 坐标.items())]
    assert len(命中) == 1, '中心在 %r 的面应恰好一个，实际 %d 个' % (坐标, len(命中))
    return 命中[0]


def 前三步(几何文件, 名, 目录, 改=None):
    '''跑第 2 步和第 3 步，返回第 2 步的结果。改 是 {面中心坐标: 类型}。'''
    标注 = None
    if 改:
        几何 = 识别几何(几何文件[名])
        标注 = {'面': {找面(几何, **{轴: 值}): 类 for (轴, 值), 类 in 改.items()}}
    结果 = generate_boundary(几何文件[名], 多孔工况() if 名 == 'C' else 基本工况(), str(目录), 标注)
    validate_boundary(结果)
    return 结果


def 划分(几何文件, 名, 目录, 改=None, 网格尺寸=None):
    边界 = 前三步(几何文件, 名, 目录, 改)
    return 边界, generate_mesh(几何文件[名], 边界, str(目录), 网格尺寸)


def 目录快照(目录):
    return {n: open(os.path.join(目录, n), 'rb').read() for n in sorted(os.listdir(目录))}


def 面积(统计):
    return {类: v['面积'] for 类, v in 统计['边界'].items()}


def 读回(路径):
    '''重新读入 mesh.msh，量出每个体的四面体数和体积、每个边界分组的面积。'''
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.open(str(路径))
        号, 坐标, _ = gmsh.model.mesh.getNodes()
        点 = {int(n): tuple(float(c) for c in 坐标[3 * i:3 * i + 3]) for i, n in enumerate(号)}
        体 = {}
        for _, t in gmsh.model.getEntities(3):
            类型, _, 节点 = gmsh.model.mesh.getElements(3, t)
            assert [int(v) for v in 类型] == [4], '体 %d 里有非四面体单元' % t
            号4 = [int(n) for n in 节点[0]]
            四面体 = [号4[i:i + 4] for i in range(0, len(号4), 4)]
            体[int(t)] = {'单元数': len(四面体), '节点': set(号4),
                          '体积': math.fsum(_四面体体积(*(点[n] for n in q)) for q in 四面体)}
        边界 = {}
        for _, g in gmsh.model.getPhysicalGroups(2):
            合计 = []
            for t in gmsh.model.getEntitiesForPhysicalGroup(2, g):
                类型, _, 节点 = gmsh.model.mesh.getElements(2, t)
                assert [int(v) for v in 类型] == [2], '面 %d 里有非三角形单元' % t
                号3 = [int(n) for n in 节点[0]]
                合计 += [_三角形面积(*(点[n] for n in 号3[i:i + 3])) for i in range(0, len(号3), 3)]
            边界[gmsh.model.getPhysicalName(2, g)] = math.fsum(合计)
        体分组 = {gmsh.model.getPhysicalName(3, g):
                  sorted(int(t) for t in gmsh.model.getEntitiesForPhysicalGroup(3, g))
                  for _, g in gmsh.model.getPhysicalGroups(3)}
        return {'点': 点, '体': 体, '边界': 边界, '体分组': 体分组}
    finally:
        gmsh.finalize()


def 中间体(边界):
    return [v['编号'] for v in 边界['几何']['体'] if abs(v['中心'][2] - 1.5) <= 绝对容差]


# ---------------------------------------------------------------- 3.1 通过

def test_T4_01(几何文件, tmp_path):
    _, 统计 = 划分(几何文件, 'A', tmp_path)
    assert 统计['网格路径'] == os.path.join(str(tmp_path), 网格文件名)
    assert os.path.getsize(统计['网格路径']) > 0
    assert 统计['单元数'] > 0
    assert 统计['网格体积'] == pytest.approx(6.0, rel=平面容差)
    assert sum(v['单元数'] for v in 读回(统计['网格路径'])['体'].values()) == 统计['单元数']


def test_T4_02(几何文件, tmp_path):
    边界, 统计 = 划分(几何文件, 'A', tmp_path)
    assert 面积(统计) == pytest.approx({'inlet': 2, 'outlet': 2, 'symmetry': 6, 'walls': 12}, rel=平面容差)
    assert sum(面积(统计).values()) == pytest.approx(22.0, rel=平面容差)
    全部面 = sorted(t for v in 统计['边界'].values() for t in v['面'])
    assert 全部面 == sorted(边界['面'])            # 每个面恰好属于一个分组


def test_T4_03(几何文件, tmp_path):
    _, 统计 = 划分(几何文件, 'C', tmp_path)
    assert len(统计['体']) == 3
    for v in 统计['体'].values():
        assert v['网格体积'] == pytest.approx(2.0, rel=平面容差)
    assert 统计['网格体积'] == pytest.approx(6.0, rel=平面容差)
    assert 面积(统计) == pytest.approx({'inlet': 2, 'outlet': 2, 'symmetry': 6, 'walls': 12}, rel=平面容差)


def test_T4_04(几何文件, tmp_path):
    边界, 统计 = 划分(几何文件, 'C', tmp_path)
    网格 = 读回(统计['网格路径'])
    点 = 网格['点']
    assert len({tuple(round(c, 9) for c in p) for p in 点.values()}) == len(点)   # 没有重复节点
    按高度 = sorted(网格['体'].values(), key=lambda v: sum(点[n][2] for n in v['节点']) / len(v['节点']))
    for 下, 上, z in ((按高度[0], 按高度[1], 1.0), (按高度[1], 按高度[2], 2.0)):
        下面 = {n for n in 下['节点'] if abs(点[n][2] - z) <= 1e-9}
        上面 = {n for n in 上['节点'] if abs(点[n][2] - z) <= 1e-9}
        assert 下面 and 下面 == 上面
    内部面 = {t for t, 类 in 边界['面'].items() if 类 == 'interface'}
    assert len(内部面) == 2
    assert not 内部面 & {t for v in 统计['边界'].values() for t in v['面']}
    assert sum(网格['边界'].values()) == pytest.approx(22.0, rel=平面容差)   # 内部面没写进文件


def test_T4_05(几何文件, tmp_path):
    边界, 统计 = 划分(几何文件, 'C', tmp_path)
    中 = 中间体(边界)
    assert len(中) == 1 and 统计['多孔体'] == 中
    assert 统计['体分组']['porous_zone'] == 中
    assert 统计['体'][中[0]]['网格体积'] == pytest.approx(2.0, rel=平面容差)
    assert sorted(统计['体分组']['fluid_zone']) == sorted(set(统计['体']) - set(中))


def test_T4_06(几何文件, tmp_path):
    边界, 统计 = 划分(几何文件, 'E', tmp_path, 弯头标注)
    几何 = 边界['几何']
    assert 统计['网格体积'] == pytest.approx(5.0, rel=平面容差)
    assert 面积(统计) == pytest.approx({'inlet': 1, 'outlet': 1, 'walls': 20}, rel=平面容差)
    assert 统计['边界']['outlet']['面'] == [找面(几何, x=4)]
    assert 找面(几何, z=0) in 统计['边界']['walls']['面']


def test_T4_07(几何文件, tmp_path):
    _, 统计 = 划分(几何文件, 'G', tmp_path, 圆柱标注, 0.05)
    assert 统计['网格体积'] == pytest.approx(math.pi / 2, rel=曲面容差)
    assert 统计['相对偏差'] <= 曲面容差
    for 类 in ('inlet', 'outlet'):
        assert 统计['边界'][类]['面积'] == pytest.approx(math.pi / 4, rel=曲面容差)


def test_T4_08(几何文件, tmp_path):
    _, 粗 = 划分(几何文件, 'A', tmp_path / 'coarse', 网格尺寸=0.5)
    _, 细 = 划分(几何文件, 'A', tmp_path / 'fine', 网格尺寸=0.25)
    for 统计 in (粗, 细):
        assert 统计['网格体积'] == pytest.approx(6.0, rel=平面容差)
    assert 细['单元数'] > 粗['单元数']
    assert (粗['网格尺寸'], 细['网格尺寸']) == (0.5, 0.25)


# ---------------------------------------------------------------- 3.2 不通过

def 圆柱粗网格报错(几何文件, 目录):
    边界 = 前三步(几何文件, 'G', 目录, 圆柱标注)
    之前 = 目录快照(str(目录))
    with pytest.raises(网格错误) as e:
        generate_mesh(几何文件['G'], 边界, str(目录), 0.5)
    assert 目录快照(str(目录)) == 之前
    return str(e.value)


def test_T4_09(几何文件, tmp_path):
    信息 = 圆柱粗网格报错(几何文件, tmp_path)
    for 词 in ('体积不符', '网格体积', '几何体积', '相对偏差'):
        assert 词 in 信息
    assert not os.path.exists(tmp_path / 网格文件名)


@pytest.mark.parametrize('尺寸', [0, -0.1, '0.5', True, float('nan'), float('inf')])
def test_T4_10(几何文件, tmp_path, 尺寸):
    边界 = 前三步(几何文件, 'A', tmp_path)
    之前 = 目录快照(str(tmp_path))
    with pytest.raises(网格错误) as e:
        generate_mesh(几何文件['A'], 边界, str(tmp_path), 尺寸)
    assert '网格尺寸' in str(e.value)
    assert 目录快照(str(tmp_path)) == 之前


def test_T4_11(几何文件, tmp_path):
    边界 = 前三步(几何文件, 'A', tmp_path)
    with pytest.raises(几何错误) as e:
        generate_mesh(str(tmp_path / 'missing.step'), 边界, str(tmp_path))
    assert 'missing.step' in str(e.value)
    坏文件, 没有面 = tmp_path / 'bad.json', tmp_path / 'noface.json'
    坏文件.write_text('{"面": {', encoding='utf-8')
    没有面.write_text('{"多孔体": []}', encoding='utf-8')
    for 输入, 词 in ((tmp_path / 'missing.json', 'missing.json'), (坏文件, 'JSON'),
                     (没有面, '面'), ({'多孔体': []}, '面')):
        with pytest.raises(标注错误) as e:
            generate_mesh(几何文件['A'], 输入, str(tmp_path / 'out'))
        assert 词 in str(e.value)
    assert not os.path.exists(tmp_path / 'out')
    assert not os.path.exists(tmp_path / 网格文件名)


def test_T4_12(几何文件, tmp_path):
    A = 前三步(几何文件, 'A', tmp_path / 'A')
    with pytest.raises(标注错误) as e:
        generate_mesh(几何文件['C'], A['标注路径'], str(tmp_path / 'out'))
    assert '外表面没有类型' in str(e.value)
    assert not os.path.exists(tmp_path / 'out')


def test_T4_13(几何文件, tmp_path):
    C = 前三步(几何文件, 'C', tmp_path / 'C')
    with pytest.raises(标注错误) as e:
        generate_mesh(几何文件['A'], C['标注路径'], str(tmp_path / 'out'))
    assert '不存在的面编号' in str(e.value)
    assert not os.path.exists(tmp_path / 'out')


# ---------------------------------------------------------------- 3.3 输入与产物

def test_T4_14(几何文件, tmp_path):
    _, 统计 = 划分(几何文件, 'C', tmp_path)
    网格 = 读回(统计['网格路径'])
    assert set(网格['体']) == set(统计['体'])
    for t, v in 统计['体'].items():
        assert 网格['体'][t]['单元数'] == v['单元数']
        assert 网格['体'][t]['体积'] == pytest.approx(v['网格体积'], rel=平面容差)
    assert 网格['边界'] == pytest.approx(面积(统计), rel=平面容差)
    assert 网格['体分组'] == {名: sorted(标签) for 名, 标签 in 统计['体分组'].items()}
    assert len(网格['点']) == 统计['节点数']


def test_T4_15(几何文件, tmp_path):
    边界 = 前三步(几何文件, 'C', tmp_path)
    甲 = generate_mesh(几何文件['C'], 边界, str(tmp_path / 'a'))
    乙 = generate_mesh(几何文件['C'], 边界['标注路径'], str(tmp_path / 'b'))
    for 键 in ('单元数', '节点数', '体分组', '多孔体'):
        assert 甲[键] == 乙[键]
    assert 面积(甲) == pytest.approx(面积(乙), rel=平面容差)
    assert {t: v['网格体积'] for t, v in 甲['体'].items()} == pytest.approx(
        {t: v['网格体积'] for t, v in 乙['体'].items()}, rel=平面容差)


def test_T4_16(几何文件, tmp_path):
    for 名, 改, 尺寸, 应通过 in (('A', None, None, True), ('G', 圆柱标注, 0.5, False)):
        目录 = str(tmp_path / 名)
        边界 = 前三步(几何文件, 名, 目录, 改)
        原边界, 之前 = copy.deepcopy(边界), 目录快照(目录)
        原几何 = open(几何文件[名], 'rb').read()
        for 输入 in (边界, 边界['标注路径']):
            if 应通过:
                generate_mesh(几何文件[名], 输入, 目录, 尺寸)
            else:
                with pytest.raises(网格错误):
                    generate_mesh(几何文件[名], 输入, 目录, 尺寸)
        之后 = 目录快照(目录)
        assert 之后[标注文件名] == 之前[标注文件名]
        assert set(之后) - set(之前) == ({网格文件名} if 应通过 else set())
        assert 边界 == 原边界
        assert open(几何文件[名], 'rb').read() == 原几何


def test_T4_17(几何文件, tmp_path):
    旧 = tmp_path / 网格文件名
    旧.write_bytes(b'old mesh')
    圆柱粗网格报错(几何文件, tmp_path)
    assert 旧.read_bytes() == b'old mesh'
    边界 = 前三步(几何文件, 'A', tmp_path)
    generate_mesh(几何文件['A'], 边界, str(tmp_path))
    assert 旧.read_bytes() != b'old mesh'


# ---------------------------------------------------------------- 3.4 可重复性

def test_T4_18(几何文件, tmp_path):
    边界 = 前三步(几何文件, 'C', tmp_path)
    甲 = generate_mesh(几何文件['C'], 边界, str(tmp_path), 0.3)
    乙 = generate_mesh(几何文件['C'], 边界, str(tmp_path), 0.3)
    assert (甲['单元数'], 甲['节点数']) == (乙['单元数'], 乙['节点数'])
    assert 甲['网格体积'] == pytest.approx(乙['网格体积'], rel=平面容差)
    assert 面积(甲) == pytest.approx(面积(乙), rel=平面容差)


def test_T4_19(几何文件, tmp_path):
    圆柱粗网格报错(几何文件, tmp_path / 'bad')
    边界, 统计 = 划分(几何文件, 'A', tmp_path / 'ok')
    assert 统计['网格体积'] == pytest.approx(6.0, rel=平面容差)
    assert 面积(统计) == pytest.approx({'inlet': 2, 'outlet': 2, 'symmetry': 6, 'walls': 12}, rel=平面容差)
