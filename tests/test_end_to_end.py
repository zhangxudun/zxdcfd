# -*- coding: utf-8 -*-
'''端到端测试: 第 1 步生成的 S1、S2 直接交给第 2 步，中间不经过人工。

链路: generate_step -> 模型.step -> generate_boundary -> 标注.json + 工况。
测试点 E2E-01 至 E2E-03 对应 test_E2E_01 至 test_E2E_03，每个都对 S1、S2 各跑一次。
预期值按 docs/01-01-生成三维结构.md 结构卡片里的公式由参数算出。面和体按位置去找，不写死编号。
运行: pip install gmsh pytest，然后在仓库根目录执行 pytest。
'''
import collections
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_structure import generate_step, 默认参数, 结构清单  # noqa: E402
from generate_boundary import generate_boundary, 标注文件名  # noqa: E402

相对容差 = 1e-6
绝对容差 = 1e-3
全部结构 = list(结构清单)


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


预期 = {'S1': 预期_S1, 'S2': 预期_S2}


# ---------------------------------------------------------------- 工具

def 近似(实际, 预期值):
    return abs(实际 - 预期值) <= 相对容差 * abs(预期值)


def 位置一致(实际, 预期值):
    return all(abs(实际[i] - 预期值[i]) <= 绝对容差 for i in range(3))


def 按类型(结果, 类):
    return [f for f in 结果['几何']['面'] if 结果['面'][f['编号']] == 类]


def 读文件(路径):
    with open(路径, encoding='utf-8') as f:
        return json.load(f)


def 跑两步(结构号, 参数, 目录):
    '''第 1 步生成 STEP，第 2 步直接读这个 STEP。返回 (STEP 路径, 第 2 步结果)。'''
    模型 = generate_step(结构号, 参数, 目录 / 'geo', 'model')
    assert os.path.isfile(模型) and os.path.getsize(模型) > 0
    return 模型, generate_boundary(模型, 多孔工况(), 目录 / 'bc')


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


# ---------------------------------------------------------------- 测试

@pytest.mark.parametrize('结构号', 全部结构)
def test_E2E_01(结构号, tmp_path):
    '''默认参数: 生成后不加标注，面类型、入口、出口、多孔体、对称面全部自动认对，工况通过。'''
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
    assert len(文件['面']) == 14
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
    assert abs(按类型(结果, 'inlet')[0]['中心'][2] - 8.0) <= 绝对容差
    assert 近似(卡['管束'][1], {'S1': 480.0, 'S2': 240.0}[结构号])
