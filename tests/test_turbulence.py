# -*- coding: utf-8 -*-
'''湍流模型（RANS: k-ε、realizable k-ε、k-ω SST，壁面函数）的测试。

TB-01  工况校验: 模型名只能是 层流 / k-ε / realizable k-ε / k-ω SST；选了模型时不能再给常数湍流粘度；
       湍流强度在 (0, 1)；层流时不能带别的键；不认识的键报错；没有 湍流 项时仍是层流。
TB-02  湍流参数: 入口当量直径 4A/P、雷诺数、k / ε / ω 的入口公式；默认混合长度 = 0.07 × 当量直径；流态提示。
TB-03  写算例: k-ω SST 有 0/k、0/omega、0/nut，壁面函数，RAS 设置，alphaD/alphaDt，wallDist，yPlus；
       k-ε 有 0/epsilon；滑移壁面不用壁面函数；定热流按 exprMixed 逐面算梯度；层流的算例文件与以前完全相同。
TB-04  第 7 步: 读 yPlus.dat，按 y+ 给近壁说明；层流时 结果.json 的格式不变。
TB-05  网页工具: 默认工况用 k-ω SST；旧版 工况.json（常数湍流粘度、没有 湍流 项）读进来仍是层流；
       整理工况去掉用不到的键。
TB-06  真实 OpenFOAM: 空气在 0.1 m 方管里 10 m/s（Re ≈ 6.4e4），k-ω SST + 壁面函数，外壁定热流:
       充分发展段的压力梯度与光滑管 Haaland 公式一致（15%），换热量与 q·A 一致（3%），壁面 y+ 在 11–300，
       第 7 步全部通过。
TB-07  真实 OpenFOAM: 水流过 +x 方向的多孔床（与 GF-10 同一几何和阻力），k-ε；第 7 步含一维对照全部通过。
'''
import copy
import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_boundary import generate_boundary, 校验工况, 工况错误, 湍流模型名  # noqa: E402
from generate_case import generate_case, 流体参数, 湍流参数, _填, _截面监测项, _向量  # noqa: E402
from generate_mesh import generate_mesh  # noqa: E402
import check_result as 检查模块  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402
import ui_tools  # noqa: E402
from test_generic_fluid import 写STEP, 水工况, 床工况, 能求解, report, x床  # noqa: E402,F401

空气 = {'工作压力': 101325.0, '摩尔质量': 0.028965, '比热': 1006.0, '分子普朗特数': 0.71}


def 管工况(湍流=None, **改):
    值 = {'入口': {'速度': 10.0, '温度': 300.0}, '出口': {'表压': 0.0}, '物性': dict(空气),
          '流向': [1, 0, 0], '壁面': '无滑移', '湍流': 湍流 or {'模型': 'k-ω SST'}}
    值.update(改)
    return 值


# ---------------------------------------------------------------- 不需要 OpenFOAM

def test_TB_01_校验():
    结果 = 校验工况(管工况(), False)
    assert 结果['湍流'] == {'模型': 'k-ω SST', '湍流强度': 0.05}
    assert 校验工况(结果, False) == 结果
    assert 湍流模型名(结果) == 'kOmegaSST'
    assert 湍流模型名(校验工况(管工况({'模型': 'k-ε'}), False)) == 'kEpsilon'
    assert 湍流模型名(校验工况(管工况({'模型': 'realizable k-ε'}), False)) == 'realizableKE'
    assert 湍流模型名(校验工况(管工况({'模型': '层流'}), False)) is None
    无 = 管工况()
    del 无['湍流']
    assert 湍流模型名(校验工况(无, False)) is None
    for 改, 关键字 in (({'湍流': {'模型': 'LES'}}, '湍流.模型'),
                       ({'湍流': {'模型': 'k-ω SST', '湍流强度': 1.0}}, '小于 1'),
                       ({'湍流': {'模型': 'k-ω SST', '湍流强度': 0}}, '湍流强度'),
                       ({'湍流': {'模型': 'k-ω SST', '混合长度': -1}}, '混合长度'),
                       ({'湍流': {'模型': 'k-ω SST', '长度': 1}}, '不认识'),
                       ({'湍流': {'模型': '层流', '湍流强度': 0.05}}, '不用填'),
                       ({'湍流': 'SST'}, '字典'),
                       ({'物性': dict(空气, 湍流粘度=1e-3)}, '湍流粘度')):
        with pytest.raises(工况错误, match=关键字):
            校验工况(管工况(**改), False)
    # 层流加常数湍流粘度仍然可以
    校验工况(管工况({'模型': '层流'}, 物性=dict(空气, 湍流粘度=1e-3)), False)


def test_TB_02_湍流参数():
    工况 = 校验工况(管工况(), False)
    流体 = 流体参数(工况)
    p = 湍流参数(工况, 流体, 0.1 * 0.2, 0.6)               # 0.1 m × 0.2 m 的矩形入口
    assert p['入口当量直径'] == pytest.approx(4 * 0.02 / 0.6)
    ν = 流体['分子粘度'] / 流体['密度']
    assert p['入口雷诺数'] == pytest.approx(10 * p['入口当量直径'] / ν)
    k = 1.5 * (10 * 0.05) ** 2
    l = 0.07 * p['入口当量直径']
    assert p['k'] == pytest.approx(k) and p['混合长度'] == pytest.approx(l)
    assert p['epsilon'] == pytest.approx(0.09 ** 0.75 * k ** 1.5 / l)
    assert p['omega'] == pytest.approx(math.sqrt(k) / (0.09 ** 0.25 * l))
    assert p['流态提示'] is None and p['混合长度来源'] == '0.07×入口当量直径'
    自给 = 校验工况(管工况({'模型': 'k-ε', '湍流强度': 0.1, '混合长度': 0.02}), False)
    q = 湍流参数(自给, 流体, 0.02, 0.6)
    assert q['混合长度'] == 0.02 and q['k'] == pytest.approx(1.5) and q['混合长度来源'] == '工况'
    慢 = 校验工况(管工况(入口={'速度': 0.1, '温度': 300.0}), False)
    assert '小于 2300' in 湍流参数(慢, 流体参数(慢), 0.02, 0.6)['流态提示']
    层 = 校验工况(管工况({'模型': '层流'}), False)
    assert '大于 4000' in 湍流参数(层, 流体参数(层), 0.02, 0.6)['流态提示']


@pytest.fixture(scope='module')
def 管(tmp_path_factory):
    目录 = tmp_path_factory.mktemp('duct')
    step = 写STEP(目录 / 'duct.step', lambda occ: occ.addBox(0, 0, 0, 2, 0.1, 0.1))
    return step, 目录


def _写(管, 工况, 名, h=0.05, 步数=20):
    step, 目录 = 管
    边界 = generate_boundary(step, 工况, 目录 / (名 + '-b'))
    validate_boundary(边界)
    网格 = generate_mesh(step, 边界, 目录 / (名 + '-m'), h)
    信息 = generate_case(step, 边界, 网格, 目录 / (名 + '-c'), 步数=步数)
    return 信息, 目录 / (名 + '-c'), 网格


def 读(路径):
    return open(路径, encoding='utf-8').read()


def test_TB_03_写算例(管):
    q = 1000.0
    信息, c, _ = _写(管, 管工况(壁面换热={'方式': '定热流', '热流密度': q}), 'sst')
    assert 信息['湍流']['OpenFOAM模型'] == 'kOmegaSST'
    assert 信息['湍流']['入口当量直径'] == pytest.approx(0.1)          # 方管 4A/P = 边长
    assert sorted(os.listdir(c / '0')) == ['T', 'U', 'k', 'nut', 'omega', 'p']
    assert 'RASModel        kOmegaSST;' in 读(c / 'constant/turbulenceProperties')
    assert 'type kqRWallFunction' in 读(c / '0/k') and 'type turbulentIntensityKineticEnergyInlet' in 读(c / '0/k')
    assert 'type omegaWallFunction' in 读(c / '0/omega')
    assert 'type turbulentMixingLengthFrequencyInlet; mixingLength 0.007;' in 读(c / '0/omega')
    assert 'type nutkWallFunction' in 读(c / '0/nut')
    控制 = 读(c / 'system/controlDict')
    流体 = 流体参数(校验工况(管工况(), False))
    α = 流体['分子导热系数'] / (流体['密度'] * 流体['比热'])
    assert 'alphaD          %.12g;' % (α / (流体['分子粘度'] / 流体['密度'])) in 控制
    assert 'alphaDt         %.12g;' % (1 / 0.9) in 控制 and '        D ' not in 控制
    assert 'type            yPlus;' in 控制
    格式 = 读(c / 'system/fvSchemes')
    assert 'div(phi,omega)' in 格式 and 'wallDist' in 格式 and 'div(phi,epsilon)' not in 格式
    assert '"(U|k|epsilon|omega)" 0.7' in 读(c / 'system/fvSolution')
    T = 读(c / '0/T')
    容 = 流体['密度'] * 流体['比热']
    assert 'type exprMixed;' in T and 'gradientExpr "%.12g/(%.12g + %.12g*nut)"' % (q / 容, α, 1 / 0.9) in T

    信息, c, _ = _写(管, 管工况({'模型': 'k-ε'}, 壁面='滑移'), 'ke')
    assert sorted(os.listdir(c / '0')) == ['T', 'U', 'epsilon', 'k', 'nut', 'p']
    assert 'WallFunction' not in 读(c / '0/k') + 读(c / '0/epsilon') + 读(c / '0/nut')
    assert 'yPlus' not in 读(c / 'system/controlDict')               # 滑移壁面没有 y+

    _, c, _ = _写(管, 管工况({'模型': '层流'}), 'lam')
    assert sorted(os.listdir(c / '0')) == ['T', 'U', 'p']
    assert 读(c / 'constant/turbulenceProperties').endswith('simulationType  laminar;\n')
    assert 'wallDist' not in 读(c / 'system/fvSchemes') and 'div(phi,k)' not in 读(c / 'system/fvSchemes')
    assert '    U\n    {' in 读(c / 'system/fvSolution') and '        D   ' in 读(c / 'system/controlDict')


def test_TB_04_yPlus(tmp_path):
    assert 检查模块._湍流结果(tmp_path, {'湍流': {'模型': '层流'}}) == {}
    assert 检查模块._湍流结果(tmp_path, {}) == {}
    结果 = 检查模块._湍流结果(tmp_path, {'湍流': {'模型': 'k-ω SST'}})
    assert 结果['壁面yPlus'] is None and '没有壁面 y+' in 结果['近壁说明']
    d = tmp_path / 'postProcessing' / 'yPlus' / '0'
    d.mkdir(parents=True)
    for 平均, 关键字 in ((5.0, '小于 11'), (60.0, '11–300'), (500.0, '大于 300')):
        (d / 'yPlus.dat').write_text('# yPlus\n# Time patch min max average\n'
                                     '1 walls 0.5 9.0 3.0\n2 walls 1.0 %g %g\n' % (平均 * 2, 平均),
                                     encoding='utf-8')
        结果 = 检查模块._湍流结果(tmp_path, {'湍流': {'模型': 'k-ε'}})
        assert 结果['壁面yPlus'] == {'最小': 1.0, '最大': 平均 * 2, '平均': 平均}
        assert 关键字 in 结果['近壁说明'] and 结果['湍流模型'] == 'k-ε'


def test_TB_05_网页工具():
    默认 = ui_tools.默认工况()
    assert 默认['湍流']['模型'] == 'k-ω SST'
    校验工况(ui_tools.整理工况(默认, True), True)
    旧 = {'物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0, '湍流粘度': 1e-3}}
    新, 不认识 = ui_tools.合并工况(默认, 旧)
    assert not 不认识 and 新['湍流'] == {'模型': '层流'}
    校验工况(ui_tools.整理工况(新, True), True)
    新, 不认识 = ui_tools.合并工况(默认, {'湍流': {'模型': 'k-ε', '混合长度': 0.3}})
    assert not 不认识 and 新['湍流'] == {'模型': 'k-ε', '混合长度': 0.3}
    整 = ui_tools.整理工况(dict(默认, 湍流={'模型': '层流', '湍流强度': 0.05}), True)
    assert 整['湍流'] == {'模型': '层流'}
    整 = ui_tools.整理工况(dict(默认, 物性=dict(默认['物性'], 湍流粘度=0.0)), True)
    assert '湍流粘度' not in 整['物性']


# ---------------------------------------------------------------- 真实 OpenFOAM

def haaland(Re):
    '''光滑管 Darcy 摩擦系数（Haaland 公式）'''
    return (-1.8 * math.log10(6.9 / Re)) ** -2


截面 = {'secA': 1.0, 'secB': 1.8}      # m，在 10D 以后的充分发展段取压


def _加截面(c):
    路径 = c / 'system/controlDict'
    文本 = 读(路径)
    加 = ''.join(_填(_截面监测项, 名=名, 点=_向量((x, 0.05, 0.05)), 法向=_向量((1, 0, 0)))
                for 名, x in 截面.items())
    i = 文本.rstrip().rfind('}')
    with open(路径, 'w', encoding='utf-8', newline='\n') as f:
        f.write(文本[:i] + 加 + 文本[i:])


def test_TB_06_方管湍流(能求解, 管):
    from run_case import run_case
    q = 1000.0
    工况 = 管工况(壁面换热={'方式': '定热流', '热流密度': q})
    信息, c, 网格 = _写(管, 工况, 'of-sst', h=0.0125, 步数=1500)
    _加截面(c)
    结果 = 检查模块.check_result(run_case(信息))
    assert 结果['通过'] is True, report(结果)
    流体 = 流体参数(校验工况(工况, False))
    Re = 信息['湍流']['入口雷诺数']
    预期 = haaland(Re) / 0.1 * 0.5 * 流体['密度'] * 10.0 ** 2           # Pa/m
    pA, pB = (检查模块._读监测(str(c), 名)[-1] * 流体['密度'] for 名 in 截面)
    梯度 = (pA - pB) / (截面['secB'] - 截面['secA'])
    说明 = 'Re %.4g，压力梯度 %.4g Pa/m，Haaland %.4g Pa/m，y+ %s' % (Re, 梯度, 预期, 结果.get('壁面yPlus'))
    print(说明)
    assert 梯度 == pytest.approx(预期, rel=0.15), 说明
    壁面积 = sum(v['面积'] for k, v in 网格['边界'].items() if k == 'walls')
    assert -结果['吸热量'] == pytest.approx(q * 壁面积, rel=0.03), report(结果)
    assert 11 < 结果['壁面yPlus']['平均'] < 300, report(结果)
    assert 结果['湍流模型'] == 'k-ω SST'


def test_TB_07_多孔床_k_epsilon(能求解, x床):
    from run_case import run_case
    step, 目录 = x床
    工况 = 床工况(湍流={'模型': 'k-ε'})
    工况['多孔'].update({'粘性阻力系数': 1e5, '惯性阻力系数': 10.0})
    边界 = generate_boundary(step, 工况, 目录 / 'tb-b')
    网格 = generate_mesh(step, 边界, 目录 / 'tb-m', 0.04)
    信息 = generate_case(step, 边界, 网格, 目录 / 'tb-c')
    assert 信息['湍流']['OpenFOAM模型'] == 'kEpsilon'
    结果 = 检查模块.check_result(run_case(信息), 对照一维=True)
    assert 结果['通过'] is True, report(结果)
    assert 300 < 结果['出口温度'] < 350
