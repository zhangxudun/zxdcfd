# -*- coding: utf-8 -*-
'''变物性（密度、粘度、导热随温度变化）与重力（浮升力）的测试，见 docs/12 的“变物性与浮升力”。

VP-01  工况校验: 物性模型 只能是 常物性 / 变物性；常物性时不能给 重力 和 温度表；重力 是三个数，默认 [0, 0, -9.81]；
       温度表 的温度严格递增、列等长、至少一列物性、不认识的列报错；表的列与对应的常数（或理想气体、Sutherland、
       分子普朗特数）只能二选一；变物性不能有常数湍流粘度；不写 物性模型 时校验结果与以前完全相同。
VP-02  物性拟合: Sutherland 的 4 次多项式在拟合范围内误差 < 0.1%；κ = μ·cp/Pr；理想气体 ρ = pM/(RT)；
       温度表 5 个点时过每个点；表没有覆盖算例温度、拟合出非正值时报工况错误。
VP-03  写算例: buoyantSimpleFoam、limitTemperature（按温度表时上下限是表的范围）、U 初场取入口速度、thermophysicalProperties（incompressiblePerfectGas / icoPolynomial，
       polynomial、hPolynomial）、constant/g、p 与 p_rgh、温差和水侧的源项作用在 h 上（codeAddSupRho），
       定热流用 externalWallHeatFluxTemperature，湍流时有 0/alphat 和 Prt；监测量多出密度；常物性算例文件不变。
VP-04  一维估算: 出口温度与常物性相同（能量方程只含 G·cp）；压降按 ρ(T)、μ(T) 沿程积分，另加加速压降
       G²(1/ρ出 − 1/ρ进)；等温时与常物性的估算相同；水侧耦合的沿程积分出口温度与 ε-NTU 一致。
VP-05  第 7 步换算（造数）: 绝对压力和质量流量换成常物性的单位，压差扣掉 ρ̄·g·Δx；日志用 log.buoyantSimpleFoam；
       没有换热时 h 的残差不参与。
VP-11  网页工具: 上传带 温度表 的物性时与之冲突的常数被删掉，上传常数时删掉表里对应的列；整理工况 在常物性时去掉 重力 和 温度表，
       变物性时去掉 湍流粘度；带温度表的预设切到变物性后能通过校验。
VP-06  真实 OpenFOAM: 700 K 烟气竖直向上流过 2 m 多孔床，温差换热冷却到约 612 K（密度变 14%）:
       第 7 步含一维对照全部通过；出口密度与理想气体 pM/(RT出) 一致（0.5%）；静压差减去流动压差 = ρ̄·g·Δx。
VP-07  真实 OpenFOAM: 同一床，水侧耦合（逆流）: 第 7 步含一维对照（烟气和水出口温度、压差）全部通过。
VP-08  真实 OpenFOAM: 浮升力。300 K 空气 1 m/s 竖直向上，多孔区里定热源加热到约 548 K，没有阻力:
       进出口静压差与 ∫ρ(T)·g·dx + G²(1/ρ出 − 1/ρ进) 一致（2%），常物性会算成 ρ入·g·L，偏大约 30%。
VP-09  真实 OpenFOAM: 水（按温度表）0.05 m/s 流过多孔床，定热源从 300 K 加热到约 338 K，不算重力:
       第 7 步含一维对照全部通过，出口密度与表一致（0.1%）。液体加重力时一维对照不适用，见 docs/12 第 10 节。
VP-10  真实 OpenFOAM: 300 K 空气 10 m/s 进 0.1 m 方管，k-ω SST、外壁定热流 1000 W/m²、变物性:
       吸热量与 q·A 一致（3%），第 7 步全部通过。
'''
import copy
import json
import math
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_boundary import generate_boundary, 校验工况, 工况错误  # noqa: E402
from generate_case import generate_case, 流体参数, 出口温度估算, 多孔压降估算, 水侧估算  # noqa: E402
from generate_mesh import generate_mesh  # noqa: E402
import check_result as 检查模块  # noqa: E402
import variable_props as vp  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402
from test_generic_fluid import 写STEP, 能求解, report  # noqa: E402,F401

气 = {'工作压力': 101325.0, '摩尔质量': 0.0289, '比热': 1100.0, '分子普朗特数': 0.7}
水表 = {'比热': 4180.0, '温度表': {'温度': [283.15, 303.15, 323.15, 343.15, 363.15],
                                  '密度': [999.7, 995.6, 988.0, 977.8, 965.3],
                                  '分子粘度': [1.306e-3, 0.797e-3, 0.547e-3, 0.404e-3, 0.315e-3],
                                  '导热系数': [0.580, 0.615, 0.643, 0.663, 0.675]}}


def 床工况(换热=None, **改):
    值 = {'物性模型': '变物性', '重力': [-9.81, 0.0, 0.0], '入口': {'速度': 5.0, '温度': 700.0},
          '出口': {'表压': 0.0}, '物性': dict(气), '流向': [1, 0, 0], '壁面': '滑移',
          # 阻力以粘性项为主: 惯性项为主时，常物性的 simpleFoam 也系统性偏低 3–4%（见 test_GF_10 的注释和 docs/12），
          # 与变物性无关；粘性项 μ(T)/ρ(T) 随温度变化最大（约 T^1.7），最能检验变物性
          '多孔': {'粘性阻力系数': 1e6, '惯性阻力系数': 1.0, '横向倍数': 10.0,
                   '换热': 换热 or {'方式': '温差', 'h_av': 800.0, '壁温': 500.0}}}
    值.update(改)
    return 值


def 读(路径):
    return open(路径, encoding='utf-8').read()


@pytest.fixture(scope='module')
def 床(tmp_path_factory):
    '''沿 x 的三段方管（1 m、2 m、1 m），中间 2 m 是多孔体，截面 0.5 m × 0.5 m。'''
    目录 = tmp_path_factory.mktemp('vpbed')
    step = 写STEP(目录 / 'bed.step', lambda occ: [occ.addBox(x0, 0, 0, L, 0.5, 0.5)
                                                   for x0, L in ((0, 1), (1, 2), (3, 1))])
    return step, 目录


@pytest.fixture(scope='module')
def 管(tmp_path_factory):
    目录 = tmp_path_factory.mktemp('vpduct')
    return 写STEP(目录 / 'duct.step', lambda occ: occ.addBox(0, 0, 0, 2, 0.1, 0.1)), 目录


def _写(step, 目录, 工况, 名, h, 步数=20):
    边界 = generate_boundary(step, 工况, 目录 / (名 + '-b'))
    validate_boundary(边界)
    网格 = generate_mesh(step, 边界, 目录 / (名 + '-m'), h)
    return generate_case(step, 边界, 网格, 目录 / (名 + '-c'), 步数=步数), 目录 / (名 + '-c'), 网格


# ---------------------------------------------------------------- 不需要 OpenFOAM

def test_VP_01_校验():
    结果 = 校验工况(床工况(), True)
    assert 结果['物性模型'] == '变物性' and 结果['重力'] == [-9.81, 0.0, 0.0]
    assert 校验工况(结果, True) == 结果
    默认 = 床工况()
    del 默认['重力']
    assert 校验工况(默认, True)['重力'] == [0.0, 0.0, -9.81]
    常 = 床工况()
    del 常['物性模型'], 常['重力']
    assert '物性模型' not in 校验工况(常, True) and '重力' not in 校验工况(常, True)
    表 = 校验工况(床工况(物性=copy.deepcopy(水表)), True)
    assert 表['物性']['温度表']['密度'][0] == 999.7
    # 只给密度列，粘度和导热用常数
    校验工况(床工况(物性={'比热': 4180.0, '分子粘度': 1e-3, '导热系数': 0.6,
                          '温度表': {'温度': [280, 360], '密度': [1000, 965]}}), True)
    for 改, 关键字 in (
            ({'物性模型': '随温度'}, '物性模型'),
            ({'物性模型': '常物性'}, '重力'),
            ({'重力': [0, -9.81]}, '重力'),
            ({'重力': [0, 'a', 0]}, '重力'),
            ({'物性': dict(气, 湍流粘度=1e-3)}, '湍流粘度'),
            ({'物性': dict(水表, 温度表={'温度': [300, 290], '密度': [996, 998]})}, '递增'),
            ({'物性': dict(水表, 温度表={'温度': [300], '密度': [996]})}, '至少两个'),
            ({'物性': dict(水表, 温度表={'温度': [300, 350], '密度': [996]})}, '一样长'),
            ({'物性': dict(水表, 温度表={'温度': [300, 350]})}, '至少要有一列'),
            ({'物性': dict(水表, 温度表={'温度': [300, 350], '比热': [4180, 4190]})}, '不认识'),
            ({'物性': dict(水表, 温度表={'温度': [300, 350], '密度': [996, -1]})}, '密度'),
            ({'物性': dict(水表, 密度=998.0)}, '二选一'),
            ({'物性': dict(水表, 工作压力=1e5, 摩尔质量=0.018)}, '二选一'),
            ({'物性': dict(水表, 分子粘度=1e-3)}, '二选一'),
            ({'物性': dict(水表, 萨瑟兰常数=110.0)}, '二选一'),
            ({'物性': dict(水表, 分子普朗特数=7.0)}, '二选一'),
            ({'物性': {'比热': 4180.0, '温度表': {'温度': [280, 360], '密度': [1000, 965]}}}, '导热系数')):
        with pytest.raises(工况错误, match=关键字):
            校验工况(床工况(**改), True)
    常表 = 床工况(物性=copy.deepcopy(水表))
    del 常表['物性模型'], 常表['重力']
    with pytest.raises(工况错误, match='温度表'):
        校验工况(常表, True)


def test_VP_02_物性拟合():
    工况 = 校验工况(床工况(), True)
    模 = vp.物性(工况, [500.0, 700.0])
    低, 高 = 模['拟合']['分子粘度']['温度范围']
    assert (低, 高) == pytest.approx((400.0, 910.0))
    for i in range(51):
        T = 低 + (高 - 低) * i / 50
        assert 模['粘度'](T) == pytest.approx(vp.sutherland(T), rel=1e-3)
        assert 模['导热'](T) == pytest.approx(vp.sutherland(T) * 1100.0 / 0.7, rel=1e-3)
        assert 模['密度'](T) == pytest.approx(101325.0 * 0.0289 / (8.314462 * T), rel=1e-12)
    assert 模['拟合']['分子粘度']['最大相对误差'] < 1e-3
    # 入口物性与常物性的写法一致（Sutherland 精确值在入口温度）
    常 = copy.deepcopy(工况)
    del 常['物性模型'], 常['重力']
    assert 流体参数(工况)['密度'] == pytest.approx(流体参数(常)['密度'], rel=1e-12)
    assert 流体参数(工况)['分子粘度'] == pytest.approx(流体参数(常)['分子粘度'], rel=1e-4)
    # 温度表: 5 个点用 4 次多项式，过每个点
    表工况 = 校验工况(床工况(物性=copy.deepcopy(水表), 入口={'速度': 0.05, '温度': 300.0}), True)
    模 = vp.物性(表工况, [300.0, 340.0])
    for i, T in enumerate(水表['温度表']['温度']):
        for 列, f in (('密度', 模['密度']), ('分子粘度', 模['粘度']), ('导热系数', 模['导热'])):
            assert f(T) == pytest.approx(水表['温度表'][列][i], rel=1e-9)
    assert 'icoPolynomial' == 模['OpenFOAM']['状态方程'][0]
    assert vp.温度限(模) == (283.15, 363.15)                      # 求解中温度限在表的范围里
    assert vp.温度限(vp.物性(工况, [500.0, 700.0])) == (250.0, 1400.0)
    with pytest.raises(工况错误, match='没有覆盖'):
        vp.物性(表工况, [300.0, 370.0])
    with pytest.raises(工况错误, match='没有覆盖'):
        vp.物性(表工况, [280.0, 340.0])
    坏 = 校验工况(床工况(物性={'比热': 4180.0, '分子粘度': 1e-3, '导热系数': 0.6,
                                  '温度表': {'温度': [300, 310, 320, 330, 340, 350, 360],
                                             '密度': [1000, 500, 1000, 500, 1000, 500, 1000]}},
                         入口={'速度': 0.05, '温度': 330.0}), True)
    with pytest.raises(工况错误, match='温度表'):
        vp.物性(坏, [330.0])


def test_VP_03_写算例(床, 管):
    step, 目录 = 床
    信息, c, _ = _写(step, 目录, 床工况(), 'w', 0.15)
    assert 信息['物性模型'] == '变物性' and 信息['求解器'] == 'buoyantSimpleFoam' and 信息['时间步'] == 1.0
    assert 信息['参考压力'] == 101325.0 and 信息['重力'] == [-9.81, 0.0, 0.0]
    assert 信息['进出口中心']['inlet'] == pytest.approx([0, 0.25, 0.25])
    assert 信息['进出口中心']['outlet'] == pytest.approx([4, 0.25, 0.25])
    assert {'inletRho', 'outletRho', 'porousInRho', 'porousOutRho'} <= set(信息['监测量'])
    assert 'application     buoyantSimpleFoam;' in 读(c / 'system/controlDict')
    热 = 读(c / 'constant/thermophysicalProperties')
    for 词 in ('heRhoThermo', 'transport       polynomial;', 'thermo          hPolynomial;',
               'equationOfState incompressiblePerfectGas;', 'pRef 101325;', 'molWeight       28.9;',
               'CpCoeffs<8>     (1100 0 0 0 0 0 0 0);', 'energy          sensibleEnthalpy;'):
        assert 词 in 热, 词
    assert 'value           (-9.81 0 0);' in 读(c / 'constant/g')
    assert sorted(os.listdir(c / '0')) == ['T', 'U', 'p', 'p_rgh']
    assert 'type prghPressure; p uniform 101325;' in 读(c / '0/p_rgh')
    assert 'fixedFluxPressure' in 读(c / '0/p_rgh')
    源 = 读(c / 'system/fvOptions')
    assert 'explicitPorositySource' in 源 and 'fields          (h);' in 源 and 'codeAddSupRho' in 源
    assert 'lookupObject<volScalarField>("T")' in 源
    assert 'buoyantSimpleFoam > log.buoyantSimpleFoam' in 读(c / 'Allrun') and 'setpriv' in 读(c / 'Allrun')
    assert 'h               1e-8;' in 读(c / 'system/fvSolution')       # h 收敛慢，1e-6 会过早停（出口温度还在变）
    assert 'type            limitTemperature;' in 源 and 'internalField   uniform (5 0 0);' in 读(c / '0/U')
    assert 'default         leastSquares;' in 读(c / 'system/fvSchemes')  # Gauss linear 在四面体上多孔压降偏低约 5.6%
    # 水侧耦合
    _, c2, _ = _写(step, 目录, 床工况({'方式': '水侧耦合', '水入口温度': 450.0, '水质量流量': 1.0, '水比热': 4300.0,
                                     'h_av': 800.0}), 'ws', 0.15)
    源 = 读(c2 / 'system/fvOptions')
    assert 'fields          (h);' in 源 and 'codeAddSupRho' in 源 and 'Cp*(Tw - T[cs[i]])' in 源
    # 定热源: 直接是 W/m³，不编译
    _, c3, _ = _写(step, 目录, 床工况({'方式': '定热源', '体积热源': 1e5}), 'q', 0.15)
    assert 'h           (100000 0);' in 读(c3 / 'system/fvOptions') and 'setpriv' not in 读(c3 / 'Allrun')
    # 液体按温度表
    信息4, c4, _ = _写(step, 目录, 床工况({'方式': '定热源', '体积热源': 4e6}, 物性=copy.deepcopy(水表),
                                          入口={'速度': 0.05, '温度': 300.0}), 'tab', 0.15)
    assert 'equationOfState icoPolynomial;' in 读(c4 / 'constant/thermophysicalProperties')
    assert 'rhoCoeffs<8>' in 读(c4 / 'constant/thermophysicalProperties')
    assert 信息4['物性拟合']['密度']['来源'] == '温度表'
    # 湍流 + 定热流
    管step, 管目录 = 管
    _, c5, _ = _写(管step, 管目录, {'物性模型': '变物性', '入口': {'速度': 10.0, '温度': 300.0}, '出口': {'表压': 0.0},
                                    '物性': dict(气), '流向': [1, 0, 0], '壁面': '无滑移',
                                    '湍流': {'模型': 'k-ω SST'}, '壁面换热': {'方式': '定热流', '热流密度': 1000.0}},
                   'sst', 0.05)
    assert sorted(os.listdir(c5 / '0')) == ['T', 'U', 'alphat', 'k', 'nut', 'omega', 'p', 'p_rgh']
    assert 'compressible::alphatWallFunction; Prt 0.9;' in 读(c5 / '0/alphat')
    assert 'kOmegaSSTCoeffs' in 读(c5 / 'constant/turbulenceProperties')
    assert 'externalWallHeatFluxTemperature; mode flux; q uniform 1000;' in 读(c5 / '0/T')
    assert os.path.isfile(c5 / 'system/snappyHexMeshDict')
    # 常物性（不写 物性模型）的算例与以前相同: simpleFoam，没有 g 和 p_rgh
    常 = 床工况()
    del 常['物性模型'], 常['重力']
    信息6, c6, _ = _写(step, 目录, 常, 'const', 0.15)
    assert 'application     simpleFoam;' in 读(c6 / 'system/controlDict')
    assert sorted(os.listdir(c6 / '0')) == ['T', 'U', 'p'] and not os.path.exists(c6 / 'constant/g')
    assert '物性模型' not in 信息6 and '求解器' not in 信息6


def test_VP_04_一维估算():
    工况 = 校验工况(床工况(), True)
    流体 = 流体参数(工况)
    G = 流体['密度'] * 5.0
    模 = vp.物性(工况, [500.0, 700.0])
    压降, 加速, T出, ρ进, ρ出 = vp.多孔压降(工况, 模, G, 2.0, 0.25)
    assert T出 == pytest.approx(出口温度估算(工况, 流体, 5.0, 2.0, 截面=0.25), rel=1e-9)
    assert ρ出 / ρ进 == pytest.approx(700.0 / T出, rel=1e-9)
    assert 加速 == pytest.approx(G * G * (1 / ρ出 - 1 / ρ进), rel=1e-12) and 加速 < 0
    常 = 多孔压降估算(工况, 流体, 5.0, 2.0)          # 常物性按入口物性
    assert 压降 < 0.95 * 常                           # 冷却后速度变慢，压降比常物性小
    # 等温时与常物性相同
    等温 = 校验工况(床工况({'方式': '定热源', '体积热源': 0.0}), True)
    p, a, _, _, _ = vp.多孔压降(等温, vp.物性(等温, [700.0]), G, 2.0, 0.25)
    assert a == pytest.approx(0, abs=1e-12)
    assert p == pytest.approx(多孔压降估算(等温, 流体参数(等温), 5.0, 2.0), rel=2e-4)
    # 水侧耦合: 沿程积分的出口温度与 ε-NTU 一致
    for 流动 in ('逆流', '顺流'):
        水 = 校验工况(床工况({'方式': '水侧耦合', '水入口温度': 450.0, '水质量流量': 1.0, '水比热': 4300.0,
                               'h_av': 800.0, '水流动': 流动}), True)
        剖面 = vp.温度剖面(水, G, 2.0, 0.25, 1100.0)
        assert 剖面[-1] == pytest.approx(水侧估算(水, 流体参数(水), 5.0, 2.0, 0.25)[0], abs=1e-6)


def _造变物性算例(root, 方程=('Ux', 'Uy', 'Uz', 'h', 'p_rgh'), 出口温度=650.0):
    from test_check_result import write_monitor
    c = Path(root) / 'vp'
    c.mkdir(parents=True)
    ρ = 0.5
    信息 = {'密度': ρ, '比热': 1100.0, '入口温度': 700.0, '出口表压': 0.0, '理论入口体积流量': 1.25,
            '多孔压降估算': 400.0, '出口温度估算': 650.0, '步数': 100, '物性模型': '变物性', '参考压力': 101325.0,
            '重力': [-9.81, 0, 0], '进出口中心': {'inlet': [0, 0.25, 0.25], 'outlet': [4, 0.25, 0.25]},
            '压力截面': {'porousInP': [1, 0, 0], 'porousOutP': [3, 0, 0]},
            '监测量': ['inletFlow', 'inletP', 'outletFlow', 'outletT', 'inletRho', 'outletRho', 'porousInP',
                       'porousOutP', 'porousInRho', 'porousOutRho']}
    (c / '算例信息.json').write_text(json.dumps(信息, ensure_ascii=False), encoding='utf-8')
    # 流动压降 400 Pa；向上流 2 m、ρ̄ 0.52: 静压再多降 0.52·9.81·2
    值 = {'inletFlow': -0.625, 'outletFlow': 0.625, 'outletT': 出口温度, 'inletRho': 0.5, 'outletRho': 0.54,
          'porousInRho': 0.5, 'porousOutRho': 0.54,
          'porousInP': 101325.0 + 410.0 + 0.52 * 9.81 * 2, 'porousOutP': 101325.0 + 10.0,
          'inletP': 101325.0 + 410.0 + 0.52 * 9.81 * 4}
    for 名, v in 值.items():
        write_monitor(c, 名, [v] * 100)
    行 = ['Starting time loop', '']
    for i in range(100):
        行 += ['Time = %d' % (i + 1), '']
        行 += ['smoothSolver:  Solving for %s, Initial residual = %s, Final residual = 1e-9, No Iterations 1'
              % (名, 0.5 if 名 == 'h' else 1e-5) for 名 in 方程]
    (c / 'log.buoyantSimpleFoam').write_text('\n'.join(行) + '\nEnd\n', encoding='utf-8')
    return c


def test_VP_05_检查结果换算(tmp_path):
    c = _造变物性算例(tmp_path / 'a', 出口温度=700.0)       # 没有换热: h 的残差 0.5 只是噪声
    结果 = 检查模块.check_result(c)
    assert 结果['入口流量'] == pytest.approx(1.25) and 结果['入口质量流量'] == pytest.approx(0.625)
    assert 结果['多孔区压差'] == pytest.approx(400.0) and 结果['多孔区平均密度'] == pytest.approx(0.52)
    assert 结果['多孔区静压差'] == pytest.approx(400.0 + 0.52 * 9.81 * 2)
    assert 结果['进出口压差'] == pytest.approx(410.0) and 结果['进出口平均密度'] == pytest.approx(0.52)
    assert 结果['出口密度'] == pytest.approx(0.54)
    assert 结果['通过'] is True, report(结果)
    assert '没有换热，h 不参与判断' in [项 for 项 in 结果['检查'] if 项['项目'] == '残差收敛'][0]['说明']
    有热 = _造变物性算例(tmp_path / 'b', 出口温度=650.0)
    结果 = 检查模块.check_result(有热, 对照一维=True)
    残差 = [项 for 项 in 结果['检查'] if 项['项目'] == '残差收敛'][0]
    assert 残差['通过'] is False and 'h 0.5' in 残差['说明']
    assert 结果['吸热量'] == pytest.approx(0.625 * 1100.0 * 50.0)
    assert [项['通过'] for 项 in 结果['检查'] if 项['项目'] == '压差与一维解一致'] == [True]


# ---------------------------------------------------------------- 真实 OpenFOAM

def _实算(step, 目录, 工况, 名, h, 步数=1500, 对照一维=True):
    from run_case import run_case
    信息, c, 网格 = _写(step, 目录, 工况, 名, h, 步数)
    运行 = run_case(信息)
    assert 'buoyantSimpleFoam' in 运行['日志']
    return 信息, 检查模块.check_result(运行, 对照一维=对照一维), c, 网格


def test_VP_11_网页工具():
    import ui_tools
    新, 不认识 = ui_tools.合并工况(ui_tools.默认工况(), {'物性模型': '变物性', '重力': [0, 0, -9.81],
                                                       '物性': copy.deepcopy(水表)})
    assert not 不认识
    for k in ('工作压力', '摩尔质量', '分子普朗特数'):
        assert k not in 新['物性']
    整 = ui_tools.整理工况(新, True)
    assert '湍流粘度' not in 整['物性'] and 整['重力'] == [0, 0, -9.81]
    校验工况(整, True)
    常 = ui_tools.整理工况(dict(新, 物性模型='常物性'), True)
    assert '重力' not in 常 and '温度表' not in 常['物性']
    回, _ = ui_tools.合并工况(新, {'物性': {'密度': 990.0}})     # 再传常数密度: 去掉表的密度列，粘度、导热仍按表
    assert '密度' not in 回['物性']['温度表'] and '分子粘度' in 回['物性']['温度表']
    校验工况(ui_tools.整理工况(回, True), True)
    物 = ui_tools.改成温度表({'密度': 1.0, '分子粘度': 1e-5, '比热': 1000.0, '分子普朗特数': 0.7},
                            {'温度': [300.0, 400.0], '分子粘度': [1.8e-5, 2.3e-5]})
    assert '分子粘度' not in 物 and 物['密度'] == 1.0 and '分子普朗特数' in 物


def test_VP_06_烟气冷却(能求解, 床):
    step, 目录 = 床
    信息, 结果, c, _ = _实算(step, 目录, 床工况(), 'of-gas', 0.04)
    print(report(结果))
    assert 结果['通过'] is True, report(结果)
    assert 结果['出口质量流量'] == pytest.approx(结果['入口质量流量'], rel=1e-4)
    ρ理 = 101325.0 * 0.0289 / (8.314462 * 结果['出口温度'])
    assert 结果['出口密度'] == pytest.approx(ρ理, rel=5e-3), report(结果)
    assert 结果['出口密度'] / 信息['密度'] > 1.1                       # 密度确实变了（约 14%）
    assert 结果['多孔区静压差'] - 结果['多孔区压差'] == pytest.approx(
        结果['多孔区平均密度'] * 9.81 * 2.0, rel=1e-3)


def test_VP_07_水侧耦合(能求解, 床):
    step, 目录 = 床
    工况 = 床工况({'方式': '水侧耦合', '水入口温度': 450.0, '水质量流量': 1.0, '水比热': 4300.0, 'h_av': 800.0})
    信息, 结果, _, _ = _实算(step, 目录, 工况, 'of-ws', 0.04)
    print(report(结果))
    assert 结果['通过'] is True, report(结果)


def test_VP_08_浮升力(能求解, 床):
    '''竖直向上的空气在中间 2 m 里被均匀加热，没有阻力: 静压差只有气柱重量和加速。'''
    step, 目录 = 床
    Q = 1.6e5
    工况 = 床工况({'方式': '定热源', '体积热源': Q}, 入口={'速度': 1.0, '温度': 300.0})
    工况['多孔'].update({'粘性阻力系数': 0.0, '惯性阻力系数': 0.0})
    信息, 结果, _, _ = _实算(step, 目录, 工况, 'of-buoy', 0.05, 对照一维=False)
    print(report(结果))
    assert 结果['通过'] is True, report(结果)
    ρ = lambda T: 101325.0 * 0.0289 / (8.314462 * T)
    G = ρ(300.0) * 1.0
    T出 = 300.0 + Q * 2.0 / (G * 1100.0)
    n = 4000
    柱 = 0.0
    for i in range(n):                      # 0–1 m 300 K，1–3 m 线性升温，3–4 m 出口温度
        x = 4.0 * (i + 0.5) / n
        T = 300.0 if x < 1 else (T出 if x > 3 else 300.0 + (T出 - 300.0) * (x - 1) / 2)
        柱 += ρ(T) * 9.81 * 4.0 / n
    预期 = 柱 + G * G * (1 / ρ(T出) - 1 / ρ(300.0))
    assert 结果['出口温度'] == pytest.approx(T出, abs=0.03 * (T出 - 300.0)), report(结果)
    assert 结果['进出口静压差'] == pytest.approx(预期, rel=0.02), (预期, report(结果))
    assert 结果['进出口静压差'] < 0.8 * ρ(300.0) * 9.81 * 4.0       # 常物性（ρ入·g·L）会偏大 30% 左右


def test_VP_09_水按温度表(能求解, 床):
    step, 目录 = 床
    工况 = 床工况({'方式': '定热源', '体积热源': 4.0e6}, 入口={'速度': 0.05, '温度': 300.0},
                  物性=copy.deepcopy(水表), 重力=[0.0, 0.0, 0.0])
    信息, 结果, _, _ = _实算(step, 目录, 工况, 'of-water', 0.04)
    print(report(结果))
    assert 结果['通过'] is True, report(结果)
    模 = vp.物性(校验工况(工况, True), [300.0, 结果['出口温度']])
    assert 结果['出口密度'] == pytest.approx(模['密度'](结果['出口温度']), rel=1e-3), report(结果)


def test_VP_10_湍流定热流(能求解, 管):
    step, 目录 = 管
    q = 1000.0
    工况 = {'物性模型': '变物性', '入口': {'速度': 10.0, '温度': 300.0}, '出口': {'表压': 0.0}, '物性': dict(气),
            '流向': [1, 0, 0], '壁面': '无滑移', '湍流': {'模型': 'k-ω SST'},
            '壁面换热': {'方式': '定热流', '热流密度': q}}
    信息, 结果, _, 网格 = _实算(step, 目录, 工况, 'of-sst', 0.025, 步数=1000, 对照一维=False)
    print(report(结果))
    assert 结果['通过'] is True, report(结果)
    壁面积 = sum(v['面积'] for k, v in 网格['边界'].items() if k == 'walls')
    assert -结果['吸热量'] == pytest.approx(q * 壁面积, rel=0.03), report(结果)
