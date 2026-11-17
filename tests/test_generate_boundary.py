# -*- coding: utf-8 -*-
'''docs/02-拿到边界条件.md 的测试，测试点 T2-01 至 T2-16 对应 test_T2_01 至 test_T2_16。

运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
测试几何 A 至 E 在测试里用 gmsh 现建并写成 STEP。面和体的编号按中心坐标去找，不写死。
'''
import collections
import copy
import json
import os
import sys
import tempfile

import gmsh
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_boundary import (generate_boundary, 识别几何, 面类型, 标注文件名,  # noqa: E402
                               标注错误, 工况错误)
from generate_structure import _写出_单位m  # noqa: E402

绝对容差 = 1e-3
_轴 = {'x': 0, 'y': 1, 'z': 2}


# ---------------------------------------------------------------- 测试几何

def _长方体(x0):
    return lambda occ: occ.addBox(x0, 0, 0, 1, 2, 3)


def _叠放(层数):
    return lambda occ: [occ.addBox(0, 0, z, 1, 2, 1) for z in range(层数)]


def _弯头(occ):
    '''L 形截面在 y=0 平面上，沿 y 拉伸 1 m。'''
    角 = [(1, 0), (4, 0), (4, 1), (2, 1), (2, 3), (1, 3)]
    点 = [occ.addPoint(x, 0, z) for x, z in 角]
    边 = [occ.addLine(点[i], 点[(i + 1) % len(点)]) for i in range(len(点))]
    截面 = occ.addPlaneSurface([occ.addCurveLoop(边)])
    occ.extrude([(2, 截面)], 0, 1, 0)


测试几何 = {'A': _长方体(0), 'B': _长方体(1), 'C': _叠放(3), 'D': _叠放(2), 'E': _弯头}


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


def 找体(几何, z):
    命中 = [v['编号'] for v in 几何['体'] if abs(v['中心'][2] - z) <= 绝对容差]
    assert len(命中) == 1, '中心在 z=%r 的体应恰好一个，实际 %d 个' % (z, len(命中))
    return 命中[0]


def 计数(结果):
    return dict(collections.Counter(结果['面'].values()))


def 读文件(路径):
    with open(路径, encoding='utf-8') as f:
        return json.load(f)


def 目录快照(目录):
    return {n: open(os.path.join(目录, n), 'rb').read() for n in sorted(os.listdir(目录))}


def 断言报错且不留产物(tmp_path, 异常, step, 工况, 标注=None):
    '''目标原来没有文件和原来有文件两种情况都测，返回报错信息。'''
    信息 = ''
    for 有旧文件 in (False, True):
        目录 = tempfile.mkdtemp(dir=str(tmp_path))
        if 有旧文件:
            with open(os.path.join(目录, 标注文件名), 'w', encoding='utf-8') as f:
                f.write('old content')
        之前 = 目录快照(目录)
        with pytest.raises(异常) as e:
            generate_boundary(step, 工况, 目录, 标注)
        assert 目录快照(目录) == 之前
        信息 = str(e.value)
    return 信息


# ---------------------------------------------------------------- 3.1 自动识别

def test_T2_01(几何文件, tmp_path):
    结果 = generate_boundary(几何文件['A'], 基本工况(), tmp_path)
    几何, 面 = 结果['几何'], 结果['面']
    assert len(面) == 6
    assert 面[找面(几何, z=3)] == 'inlet'
    assert 面[找面(几何, z=0)] == 'outlet'
    assert 面[找面(几何, x=0)] == 'symmetry'
    assert 计数(结果) == {'inlet': 1, 'outlet': 1, 'symmetry': 1, 'walls': 3}
    assert 结果['多孔体'] == []


def test_T2_02(几何文件, tmp_path):
    结果 = generate_boundary(几何文件['B'], 基本工况(), tmp_path)
    assert 计数(结果) == {'inlet': 1, 'outlet': 1, 'walls': 4}


def test_T2_03(几何文件, tmp_path):
    结果 = generate_boundary(几何文件['C'], 多孔工况(), tmp_path)
    assert 计数(结果) == {'inlet': 1, 'outlet': 1, 'interface': 2, 'symmetry': 3, 'walls': 9}
    assert 结果['多孔体'] == [找体(结果['几何'], 1.5)]


def test_T2_04(几何文件, tmp_path):
    结果 = generate_boundary(几何文件['D'], 基本工况(), tmp_path)
    assert 计数(结果) == {'inlet': 1, 'outlet': 1, 'interface': 1, 'symmetry': 2, 'walls': 6}
    assert 结果['多孔体'] == []


def test_T2_05(几何文件, tmp_path):
    '''已知限制: 弯头的出口被猜错。规则改了要同步改这一条。'''
    结果 = generate_boundary(几何文件['E'], 基本工况(), tmp_path)
    几何, 面 = 结果['几何'], 结果['面']
    assert len(面) == 8
    assert 面[找面(几何, z=3)] == 'inlet'
    assert 面[找面(几何, z=0)] == 'outlet'
    assert 面[找面(几何, x=4)] == 'walls'


# ---------------------------------------------------------------- 3.2 标注覆盖

def test_T2_06(几何文件, tmp_path):
    几何 = 识别几何(几何文件['E'])
    底, 端, 顶 = 找面(几何, z=0), 找面(几何, x=4), 找面(几何, z=3)
    结果 = generate_boundary(几何文件['E'], 基本工况(), tmp_path,
                             {'面': {底: 'walls', 端: 'outlet'}})
    assert 结果['面'][底] == 'walls' and 结果['面'][端] == 'outlet'
    assert 结果['面'][顶] == 'inlet'
    assert 计数(结果) == {'inlet': 1, 'outlet': 1, 'walls': 6}
    自动 = {f['编号']: f['类型'] for f in 几何['面']}
    for t, 类 in 结果['面'].items():
        if t not in (底, 端):
            assert 类 == 自动[t]


def test_T2_07(几何文件, tmp_path):
    几何 = 识别几何(几何文件['D'])
    下 = 找体(几何, 0.5)
    结果 = generate_boundary(几何文件['D'], 多孔工况(), tmp_path, {'多孔体': [下]})
    assert 结果['多孔体'] == [下]
    assert 结果['面'] == {f['编号']: f['类型'] for f in 几何['面']}


def test_T2_08(几何文件, tmp_path):
    几何 = 识别几何(几何文件['C'])
    上, 中 = 找体(几何, 2.5), 找体(几何, 1.5)
    assert 几何['多孔体'] == [中]
    结果 = generate_boundary(几何文件['C'], 多孔工况(), tmp_path, {'多孔体': [上]})
    assert 结果['多孔体'] == [上]


def test_T2_09(几何文件, tmp_path):
    坏文件 = str(tmp_path / 'bad.json')
    with open(坏文件, 'w', encoding='utf-8') as f:
        f.write('{"面": {')
    用例 = [({'面': {99999: 'inlet'}}, '99999'),
            ({'多孔体': [99999]}, '99999'),
            ({'面': {找面(识别几何(几何文件['A']), z=3): '定壁温'}}, '定壁温'),
            (坏文件, 'JSON')]
    for 标注, 关键词 in 用例:
        信息 = 断言报错且不留产物(tmp_path, 标注错误, 几何文件['A'], 基本工况(), 标注)
        assert 关键词 in 信息, (标注, 信息)


def test_T2_10(几何文件, tmp_path):
    几何 = 识别几何(几何文件['E'])
    标注 = {'面': {找面(几何, z=0): 'walls', 找面(几何, x=4): 'outlet'}}
    第一次 = generate_boundary(几何文件['E'], 基本工况(), tmp_path / 'a', 标注)
    文件 = 读文件(第一次['标注路径'])
    assert all(isinstance(k, str) for k in 文件['面'])
    assert {int(k): v for k, v in 文件['面'].items()} == 第一次['面']
    assert 文件['多孔体'] == 第一次['多孔体']
    第二次 = generate_boundary(几何文件['E'], 基本工况(), tmp_path / 'b', 第一次['标注路径'])
    assert 第二次['面'] == 第一次['面']
    assert 第二次['多孔体'] == 第一次['多孔体']
    assert 读文件(第二次['标注路径']) == 文件


# ---------------------------------------------------------------- 3.3 工况字典

def test_T2_11(几何文件, tmp_path):
    工况 = 基本工况()
    结果 = generate_boundary(几何文件['A'], 工况, tmp_path)
    assert os.path.basename(结果['标注路径']) == 标注文件名
    文件 = 读文件(结果['标注路径'])
    assert len(文件['面']) == 6
    assert all(类 in 面类型 for 类 in 文件['面'].values())
    assert 结果['工况'] == 工况
    assert 结果['工况']['入口'] == {'速度': 5.0, '温度': 600.0}
    assert 结果['工况']['出口'] == {'表压': 0.0}
    assert 结果['工况']['壁面'] == '滑移'


@pytest.mark.parametrize('组, 项', [('入口', '速度'), ('入口', '温度'), ('出口', '表压')])
def test_T2_12(组, 项, 几何文件, tmp_path):
    工况 = 基本工况()
    del 工况[组][项]
    信息 = 断言报错且不留产物(tmp_path, 工况错误, 几何文件['A'], 工况)
    assert '%s.%s' % (组, 项) in 信息


def test_T2_13(几何文件, tmp_path):
    工况 = 多孔工况()
    结果 = generate_boundary(几何文件['C'], 工况, tmp_path)
    assert 结果['多孔体'] == [找体(结果['几何'], 1.5)]
    assert 结果['工况']['多孔'] == 工况['多孔']
    assert 读文件(结果['标注路径'])['多孔体'] == 结果['多孔体']


@pytest.mark.parametrize('路径, 关键词', [
    (('多孔',), '多孔'),
    (('多孔', '粘性阻力系数'), '粘性阻力系数'),
    (('多孔', '惯性阻力系数'), '惯性阻力系数'),
    (('多孔', '换热'), '换热'),
    (('多孔', '换热', 'h_av'), 'h_av'),
])
def test_T2_14(路径, 关键词, 几何文件, tmp_path):
    工况 = copy.deepcopy(多孔工况())
    上级 = 工况
    for 键 in 路径[:-1]:
        上级 = 上级[键]
    del 上级[路径[-1]]
    信息 = 断言报错且不留产物(tmp_path, 工况错误, 几何文件['C'], 工况)
    assert 关键词 in 信息


@pytest.mark.parametrize('项, 值', [
    ('速度', -1.0), ('速度', 0), ('速度', 'abc'), ('速度', None), ('速度', True),
    ('速度', float('nan')), ('速度', float('inf')),
    ('温度', 0), ('温度', -10.0),
])
def test_T2_15(项, 值, 几何文件, tmp_path):
    工况 = 基本工况()
    工况['入口'][项] = 值
    信息 = 断言报错且不留产物(tmp_path, 工况错误, 几何文件['A'], 工况)
    assert '入口.' + 项 in 信息


# ---------------------------------------------------------------- 3.4 可重复性

def test_T2_16(几何文件, tmp_path):
    第一次 = generate_boundary(几何文件['C'], 多孔工况(), tmp_path)
    内容 = open(第一次['标注路径'], 'rb').read()
    第二次 = generate_boundary(几何文件['C'], 多孔工况(), tmp_path)
    assert 第二次['面'] == 第一次['面']
    assert 第二次['多孔体'] == 第一次['多孔体']
    assert open(第二次['标注路径'], 'rb').read() == 内容
    assert os.listdir(str(tmp_path)) == [标注文件名]
