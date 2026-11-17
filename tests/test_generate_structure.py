# -*- coding: utf-8 -*-
'''docs/01-01-生成三维结构.md 的测试，规则 G01 至 G15。

测试点 T1-01-<结构号>-<规则号> 对应 test_<规则号>[<结构号>]，T1-01-ALL-G15c 对应 test_G15c。
运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
所有几何量都在读回并 fragment 之后的 STEP 上量。
'''
import os
import shutil
import sys
import tempfile

import gmsh
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_structure import generate_step, 默认参数, 参数错误, 文件名错误, 结构清单  # noqa: E402

相对容差 = 1e-6
绝对容差 = 1e-3
全部结构 = list(结构清单)


# ---------------------------------------------------------------- 结构卡片

def 过渡段体积(a1, b1, a2, b2, h):
    return h / 6 * (2 * a1 * b1 + 2 * a2 * b2 + a1 * b2 + a2 * b1)


def 卡片_S1(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    return {
        '实体': {'入口段': ((0, b / 2, h + li / 2), 2 * a * b * li),
                 '管束': ((0, b / 2, h / 2), 2 * a * b * h),
                 '出口段': ((0, b / 2, -lo / 2), 2 * a * b * lo)},
        '面数': 16,
        '包围盒': ((-a, a), (0, b), (-lo, h + li)),
        '入口面': (((-a, a), (0, b), (h + li, h + li)), 2 * a * b),
        '出口面': (((-a, a), (0, b), (-lo, -lo)), 2 * a * b),
        '共用面': [(((-a, a), (0, b), (0, 0)), 2 * a * b),
                   (((-a, a), (0, b), (h, h)), 2 * a * b)],
        '对称面': (0, 0.0),
    }


def 卡片_S2(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    ai, bi = p['入口半宽x'], p['入口深度y']
    ao, bo = p['出口半宽x'], p['出口深度y']
    return {
        '实体': {'入口段': ((ai / 2, b / 2, h + li / 2), 过渡段体积(a, b, ai, bi, li)),
                 '管束': ((a / 2, b / 2, h / 2), a * b * h),
                 '出口段': ((ao / 2, b / 2, -lo / 2), 过渡段体积(a, b, ao, bo, lo))},
        '面数': 16,
        '包围盒': ((0, a), (0, b), (-lo, h + li)),
        '入口面': (((0, ai), ((b - bi) / 2, (b + bi) / 2), (h + li, h + li)), ai * bi),
        '出口面': (((0, ao), ((b - bo) / 2, (b + bo) / 2), (-lo, -lo)), ao * bo),
        '共用面': [(((0, a), (0, b), (0, 0)), a * b),
                   (((0, a), (0, b), (h, h)), a * b)],
        '对称面': (3, (b + bi) / 2 * li + b * h + (b + bo) / 2 * lo),
    }


卡片 = {'S1': 卡片_S1, 'S2': 卡片_S2}

# 卡片里的默认体积，顺序: 入口段、管束、出口段
默认体积 = {'S1': (288.0, 307.2, 384.0), 'S2': (84.0, 153.6, 284 / 3)}

# 改参数用例: (改动, 预期体积, 包围盒 z 范围, 共用面的 z)
改参数用例 = {
    'S1': ({'管束高度': 5.0}, {'入口段': 288.0, '管束': 480.0, '出口段': 384.0}, (-4.0, 8.0), (0.0, 5.0)),
    'S2': ({'管束高度': 5.0}, {'入口段': 84.0, '管束': 240.0, '出口段': 284 / 3}, (-4.0, 8.0), (0.0, 5.0)),
}

# 违反约束用例: (结构号, 参数名, 值, 报错信息里要出现的词)
违反约束用例 = [
    ('S2', '入口半宽x', 7.0, '入口'),
    ('S2', '入口深度y', 9.0, '入口'),
    ('S2', '出口半宽x', 7.0, '出口'),
    ('S2', '出口深度y', 9.0, '出口'),
]

# 边界值用例: (结构号, 改动, 入口段体积, 入口面面积)
边界值用例 = [('S2', {'入口半宽x': 6.0, '入口深度y': 8.0}, 144.0, 48.0)]


# ---------------------------------------------------------------- 读回与量取

def 读回(路径):
    '''按 m 读回 STEP，做 fragment，返回纯数据: 体、面、顶点坐标、体积、面积。'''
    临时目录 = tempfile.mkdtemp(prefix='readstep_')
    try:
        副本 = os.path.join(临时目录, 'model.step')
        shutil.copyfile(os.fspath(路径), 副本)
        gmsh.initialize()
        try:
            gmsh.option.setNumber('General.Terminal', 0)
            gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
            gmsh.model.add('readback')
            gmsh.model.occ.importShapes(副本)
            gmsh.model.occ.synchronize()
            体 = gmsh.model.getEntities(3)
            模型 = {'前实体数': len(体),
                    '前总体积': sum(gmsh.model.occ.getMass(3, t) for _, t in 体)}
            if len(体) > 1:
                gmsh.model.occ.fragment(体[:1], 体[1:])
                gmsh.model.occ.synchronize()
            面 = {}
            for _, s in gmsh.model.getEntities(2):
                点 = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False, recursive=True)
                面[s] = {'面积': gmsh.model.occ.getMass(2, s),
                         '顶点': [tuple(gmsh.model.getValue(0, t, [])) for _, t in 点],
                         '所属体': []}
            模型['体'] = []
            for 序号, (_, v) in enumerate(gmsh.model.getEntities(3)):
                边界 = gmsh.model.getBoundary([(3, v)], combined=False, oriented=False, recursive=False)
                面号 = sorted(set(abs(t) for _, t in 边界))
                for s in 面号:
                    面[s]['所属体'].append(序号)
                模型['体'].append({'体积': gmsh.model.occ.getMass(3, v), '面': 面号})
            模型['面'] = 面
            return 模型
        finally:
            gmsh.finalize()
    finally:
        shutil.rmtree(临时目录, ignore_errors=True)


def 近似(实际, 预期):
    return abs(实际 - 预期) <= 相对容差 * abs(预期)


def 范围(顶点):
    return tuple((min(c[i] for c in 顶点), max(c[i] for c in 顶点)) for i in range(3))


def 范围一致(实际, 预期):
    return all(abs(实际[i][j] - 预期[i][j]) <= 绝对容差 for i in range(3) for j in range(2))


def 包围盒(模型):
    return 范围([c for f in 模型['面'].values() for c in f['顶点']])


def 找面(模型, 预期范围):
    return [f for f in 模型['面'].values() if 范围一致(范围(f['顶点']), 预期范围)]


def _减(u, v):
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def _叉(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def _点(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def 包含点(模型, 体, 点):
    '''点是否严格在实体内部。只适用于各面都是平面的凸实体。'''
    顶点 = [c for s in 体['面'] for c in 模型['面'][s]['顶点']]
    中心 = tuple(sum(c[i] for c in 顶点) / len(顶点) for i in range(3))
    for s in 体['面']:
        v = 模型['面'][s]['顶点']
        p0 = v[0]
        p1 = max(v, key=lambda c: _点(_减(c, p0), _减(c, p0)))
        法向 = max((_叉(_减(p1, p0), _减(c, p0)) for c in v), key=lambda n: _点(n, n))
        侧 = _点(法向, _减(点, p0)) * _点(法向, _减(中心, p0))
        if 侧 <= 0:
            return False
    return True


def 找实体(模型, 点):
    命中 = [b for b in 模型['体'] if 包含点(模型, b, 点)]
    assert len(命中) == 1, '参考点 %r 应恰好落在一个实体内，实际 %d 个' % (点, len(命中))
    return 命中[0]


def 生成并读回(结构号, 参数, 目录, 文件名='model'):
    路径 = generate_step(结构号, 参数, 目录, 文件名)
    assert os.path.isfile(路径) and os.path.getsize(路径) > 0
    return 读回(路径)


def 检查基本项(模型, 卡):
    '''G01、G02、G04、G05。'''
    assert len(模型['体']) == len(卡['实体'])
    for 名, (参考点, 体积) in 卡['实体'].items():
        assert 近似(找实体(模型, 参考点)['体积'], 体积), 名
    assert len(模型['面']) == 卡['面数']
    assert 范围一致(包围盒(模型), 卡['包围盒'])


def 摘要一致(甲, 乙):
    return (len(甲['体']) == len(乙['体']) and len(甲['面']) == len(乙['面'])
            and all(近似(x, y) for x, y in zip(sorted(b['体积'] for b in 甲['体']),
                                              sorted(b['体积'] for b in 乙['体'])))
            and 范围一致(包围盒(甲), 包围盒(乙)))


def 目录快照(目录):
    return {n: open(os.path.join(目录, n), 'rb').read() for n in sorted(os.listdir(目录))}


def 断言报错且不生成文件(tmp_path, 结构号, 参数, 文件名, 异常, 预置旧文件=True):
    '''目标原来没有文件和原来有文件两种情况都测，返回报错信息。'''
    信息 = ''
    for 有旧文件 in ((False, True) if 预置旧文件 else (False,)):
        目录 = tempfile.mkdtemp(dir=str(tmp_path))
        if 有旧文件:
            with open(os.path.join(目录, 文件名 + '.step'), 'w') as f:
                f.write('old content')
        之前 = 目录快照(目录)
        with pytest.raises(异常) as e:
            generate_step(结构号, 参数, 目录, 文件名)
        assert 目录快照(目录) == 之前
        信息 = str(e.value)
    return 信息


@pytest.fixture(scope='module')
def 默认模型(tmp_path_factory):
    缓存 = {}

    def 取(结构号):
        if 结构号 not in 缓存:
            目录 = tmp_path_factory.mktemp('default_' + 结构号)
            缓存[结构号] = 生成并读回(结构号, 默认参数(结构号), 目录)
        return 缓存[结构号]
    return 取


# ---------------------------------------------------------------- G01 至 G10

@pytest.mark.parametrize('结构号', 全部结构)
def test_G01(结构号, 默认模型):
    模型 = 默认模型(结构号)
    assert len(模型['体']) == len(卡片[结构号](默认参数(结构号))['实体'])


@pytest.mark.parametrize('结构号', 全部结构)
def test_G02(结构号, 默认模型):
    模型 = 默认模型(结构号)
    实体 = 卡片[结构号](默认参数(结构号))['实体']
    for (名, (参考点, 体积)), 卡片值 in zip(实体.items(), 默认体积[结构号]):
        assert 近似(体积, 卡片值), 名
        assert 近似(找实体(模型, 参考点)['体积'], 体积), 名


@pytest.mark.parametrize('结构号', 全部结构)
def test_G03(结构号, 默认模型):
    模型 = 默认模型(结构号)
    assert len(模型['体']) == 模型['前实体数']
    assert 近似(sum(b['体积'] for b in 模型['体']), 模型['前总体积'])


@pytest.mark.parametrize('结构号', 全部结构)
def test_G04(结构号, 默认模型):
    assert len(默认模型(结构号)['面']) == 卡片[结构号](默认参数(结构号))['面数']


@pytest.mark.parametrize('结构号', 全部结构)
def test_G05(结构号, 默认模型):
    assert 范围一致(包围盒(默认模型(结构号)), 卡片[结构号](默认参数(结构号))['包围盒'])


@pytest.mark.parametrize('结构号', 全部结构)
def test_G06(结构号, 默认模型):
    模型 = 默认模型(结构号)
    卡 = 卡片[结构号](默认参数(结构号))
    for 名 in ('入口面', '出口面'):
        预期范围, 面积 = 卡[名]
        命中 = 找面(模型, 预期范围)
        assert len(命中) == 1, 名
        assert 近似(命中[0]['面积'], 面积), 名
        assert len(命中[0]['所属体']) == 1, 名


@pytest.mark.parametrize('结构号', 全部结构)
def test_G07(结构号, 默认模型):
    模型 = 默认模型(结构号)
    声明 = 卡片[结构号](默认参数(结构号))['共用面']
    共用 = [f for f in 模型['面'].values() if len(f['所属体']) == 2]
    assert len(共用) == len(声明)
    for 预期范围, 面积 in 声明:
        命中 = [f for f in 共用 if 范围一致(范围(f['顶点']), 预期范围)]
        assert len(命中) == 1
        assert 近似(命中[0]['面积'], 面积)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G08(结构号, 默认模型):
    模型 = 默认模型(结构号)
    个数, 面积合计 = 卡片[结构号](默认参数(结构号))['对称面']
    在对称面上 = [f for f in 模型['面'].values()
                  if all(abs(c[0]) <= 绝对容差 for c in f['顶点'])]
    assert len(在对称面上) == 个数
    if 个数:
        assert 近似(sum(f['面积'] for f in 在对称面上), 面积合计)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G09(结构号, tmp_path):
    路径 = generate_step(结构号, 默认参数(结构号), tmp_path, 'unit')
    with open(路径, encoding='utf-8', errors='ignore') as f:
        文本 = ''.join(f.read().split()).upper()
    片段 = 文本.split('LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(')[1:]
    assert 片段, 'STEP 文件里没有找到长度单位声明'
    for 段 in 片段:
        assert 段.split(')')[0] == '$,.METRE.', '长度单位不是 m: ' + 段.split(')')[0]
    assert 范围一致(包围盒(读回(路径)), 卡片[结构号](默认参数(结构号))['包围盒'])


@pytest.mark.parametrize('结构号', 全部结构)
def test_G10(结构号, tmp_path):
    改动, 体积, z范围, 共用z = 改参数用例[结构号]
    参数 = dict(默认参数(结构号), **改动)
    卡 = 卡片[结构号](参数)
    模型 = 生成并读回(结构号, 参数, tmp_path)
    for 名, (参考点, 公式值) in 卡['实体'].items():
        assert 近似(公式值, 体积[名]), 名
        assert 近似(找实体(模型, 参考点)['体积'], 体积[名]), 名
    盒 = 包围盒(模型)
    assert 范围一致(盒, 卡['包围盒'])
    assert abs(盒[2][0] - z范围[0]) <= 绝对容差 and abs(盒[2][1] - z范围[1]) <= 绝对容差
    共用 = [f for f in 模型['面'].values() if len(f['所属体']) == 2]
    assert len(共用) == len(共用z)
    实际z = sorted(范围(f['顶点'])[2] for f in 共用)
    for (低, 高), z in zip(实际z, sorted(共用z)):
        assert abs(低 - z) <= 绝对容差 and abs(高 - z) <= 绝对容差


# ---------------------------------------------------------------- G11 至 G14

@pytest.mark.parametrize('结构号', 全部结构)
def test_G11a(结构号, tmp_path):
    for 名 in 默认参数(结构号):
        参数 = 默认参数(结构号)
        del 参数[名]
        assert 名 in 断言报错且不生成文件(tmp_path, 结构号, 参数, 'model', 参数错误)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G11b(结构号, tmp_path):
    非法值 = [0, -1.0, 'abc', '3.0', True, None, float('nan'), float('inf'), float('-inf')]
    for 名 in 默认参数(结构号):
        for 值 in 非法值:
            参数 = dict(默认参数(结构号), **{名: 值})
            assert 名 in 断言报错且不生成文件(tmp_path, 结构号, 参数, 'model', 参数错误), (名, 值)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G11c(结构号, tmp_path):
    参数 = dict(默认参数(结构号), 多余参数=1.0)
    assert '多余参数' in 断言报错且不生成文件(tmp_path, 结构号, 参数, 'model', 参数错误)


@pytest.mark.parametrize('结构号, 名, 值, 关键词', 违反约束用例)
def test_G12(结构号, 名, 值, 关键词, tmp_path):
    参数 = dict(默认参数(结构号), **{名: 值})
    信息 = 断言报错且不生成文件(tmp_path, 结构号, 参数, 'model', 参数错误)
    assert 关键词 in 信息 and 名 in 信息


@pytest.mark.parametrize('结构号, 改动, 入口段体积, 入口面面积', 边界值用例)
def test_G13(结构号, 改动, 入口段体积, 入口面面积, tmp_path):
    参数 = dict(默认参数(结构号), **改动)
    卡 = 卡片[结构号](参数)
    模型 = 生成并读回(结构号, 参数, tmp_path)
    检查基本项(模型, 卡)
    assert 近似(找实体(模型, 卡['实体']['入口段'][0])['体积'], 入口段体积)
    命中 = 找面(模型, 卡['入口面'][0])
    assert len(命中) == 1 and 近似(命中[0]['面积'], 入口面面积)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G14a(结构号, tmp_path):
    for 文件名 in ('管道', 'a b', 'a*b', 'a/b', 'a.step', ''):
        断言报错且不生成文件(tmp_path, 结构号, 默认参数(结构号), 文件名, 文件名错误, 预置旧文件=False)


@pytest.mark.parametrize('结构号', 全部结构)
def test_G14b(结构号, tmp_path):
    for 目录 in (tmp_path, tmp_path / '中文目录', tmp_path / 'new' / 'nested'):
        路径 = generate_step(结构号, 默认参数(结构号), 目录, 'test_duct-01')
        assert os.path.basename(路径) == 'test_duct-01.step'
        assert os.path.isfile(路径) and os.path.getsize(路径) > 0


# ---------------------------------------------------------------- G15

@pytest.mark.parametrize('结构号', 全部结构)
def test_G15a(结构号, tmp_path):
    卡 = 卡片[结构号](默认参数(结构号))
    第一次 = 生成并读回(结构号, 默认参数(结构号), tmp_path, 'repeat')
    第二次 = 生成并读回(结构号, 默认参数(结构号), tmp_path, 'repeat')
    检查基本项(第一次, 卡)
    检查基本项(第二次, 卡)
    assert 摘要一致(第一次, 第二次)
    assert os.listdir(str(tmp_path)) == ['repeat.step']


@pytest.mark.parametrize('结构号', 全部结构)
def test_G15b(结构号, tmp_path):
    with pytest.raises(参数错误):
        generate_step(结构号, dict(默认参数(结构号), 半宽x=-1.0), tmp_path, 'model')
    检查基本项(生成并读回(结构号, 默认参数(结构号), tmp_path), 卡片[结构号](默认参数(结构号)))


def test_G15c(tmp_path, 默认模型):
    for 序号, 结构号 in enumerate(全部结构 + 全部结构[:1]):
        模型 = 生成并读回(结构号, 默认参数(结构号), tmp_path, 'seq_%d' % 序号)
        检查基本项(模型, 卡片[结构号](默认参数(结构号)))
        assert 摘要一致(模型, 默认模型(结构号))
