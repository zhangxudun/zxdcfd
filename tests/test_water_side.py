# -*- coding: utf-8 -*-
'''多孔区“水侧耦合”换热的测试（管内工质温度由能量平衡算出，不再人工给壁温）。

WS-01  工况校验: h_av 与 管束 二选一；管内径 < 管外径；水流动 只能逆流 / 顺流；分段数 是 1–1000 的整数；
       校验两次结果不变（不把算出的 h_av 写回工况）。
WS-02  管束换热系数: 烟气侧、污垢、管壁（圆筒导热）、水侧（折算到外表面）四项热阻串联。
WS-03  一维 ε-NTU 估算: 烟气放热 = 水吸热；逆流比顺流换热多；Cr = 1 的逆流公式；NTU 很大时逆流逼近极限。
WS-04  写算例: 温度方程 fvOptions 是 scalarCodedSource，常数写对；算例信息.json 里有 水侧 一项；伪时间步受 1/k 限制。
WS-05  第 7 步（造数）: 水侧三项检查和水出口一维对照的顺序与判定；外壁非绝热时热平衡不对照；饱和温度。
WS-06  第 8 步: 认水侧耦合的七项 / 十项检查；水未沸腾不通过时回到第 2 步。
WS-07  真实 OpenFOAM: 空气横掠 +x 方向的 2 m 管束（按管束尺寸算 h_av），水逆流；第 7 步含一维对照全部通过。
WS-08  真实 OpenFOAM: 同一管束，水顺流、总 h_av 直接给；第 7 步含一维对照全部通过，水出口比逆流低。
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
from generate_boundary import generate_boundary, 校验工况, 工况错误, 管束换热系数, 体积换热系数  # noqa: E402
from generate_case import generate_case, 流体参数, 水侧估算, 伪时间步  # noqa: E402
from generate_mesh import generate_mesh  # noqa: E402
from check_result import check_result  # noqa: E402
from judge_result import judge_result, 判定错误  # noqa: E402
from validate_boundary import validate_boundary  # noqa: E402
import test_check_result as 造数  # noqa: E402
from test_generic_fluid import 写STEP, 能求解, report  # noqa: E402,F401

管束 = {'单位体积换热面积': 15.0, '烟气侧换热系数': 60.0, '管外径': 0.038, '管内径': 0.030,
        '管壁导热系数': 45.0, '水侧换热系数': 8000.0, '污垢热阻': 0.002}


def 工况(换热=None, **改):
    '''空气（理想气体）沿 +x 横掠管束，水在管内。'''
    换 = {'方式': '水侧耦合', '水入口温度': 450.0, '水质量流量': 0.05, '水比热': 4400.0, '管束': dict(管束)}
    换.update(换热 or {})
    值 = {'入口': {'速度': 2.0, '温度': 650.0}, '出口': {'表压': 0.0},
          '物性': {'比热': 1080.0, '工作压力': 101325.0, '摩尔质量': 0.02896},
          '流向': [1, 0, 0], '壁面': '滑移',
          '多孔': {'粘性阻力系数': 1e4, '惯性阻力系数': 10.0, '横向倍数': 10.0, '换热': 换}}
    值.update(改)
    return 值


def test_WS_01_校验():
    结果 = 校验工况(工况(), True)
    换 = 结果['多孔']['换热']
    assert 换['水流动'] == '逆流' and 换['分段数'] == 40 and 'h_av' not in 换
    assert 校验工况(结果, True) == 结果
    for 改, 关键字 in (({'h_av': 500.0}, '二选一'),
                       ({'管束': dict(管束, 管内径=0.04)}, '管内径'),
                       ({'水流动': '叉流'}, '水流动'),
                       ({'分段数': 0}, '分段数'),
                       ({'分段数': 2.5}, '分段数'),
                       ({'水质量流量': 0}, '水质量流量'),
                       ({'管束': {k: v for k, v in 管束.items() if k != '水侧换热系数'}}, '水侧换热系数')):
        with pytest.raises(工况错误, match=关键字):
            校验工况(工况(改), True)
    只给总系数 = 工况()
    del 只给总系数['多孔']['换热']['管束']
    只给总系数['多孔']['换热']['h_av'] = 700
    assert 体积换热系数(校验工况(只给总系数, True)['多孔']['换热']) == 700.0
    缺 = 工况()
    del 缺['多孔']['换热']['管束']
    with pytest.raises(工况错误, match='二选一'):
        校验工况(缺, True)
    assert 校验工况(工况({'水饱和温度': 580}), True)['多孔']['换热']['水饱和温度'] == 580.0


def test_WS_02_管束换热系数():
    值 = 管束换热系数(管束)
    热阻 = 值['热阻']
    assert 热阻['烟气侧'] == pytest.approx(1 / 60)
    assert 热阻['污垢'] == 0.002
    assert 热阻['管壁'] == pytest.approx(0.038 * math.log(0.038 / 0.030) / (2 * 45))
    assert 热阻['水侧'] == pytest.approx(0.038 / (0.030 * 8000))
    assert 值['总传热系数'] == pytest.approx(1 / sum(热阻.values()))
    assert 值['h_av'] == pytest.approx(15 * 值['总传热系数'])
    assert 值['总传热系数'] < 60                      # 串联热阻只会更小


def _估算(改=None):
    结果 = 校验工况(工况(改), True)
    流体 = 流体参数(结果)
    return 结果, 流体, 水侧估算(结果, 流体, 2.0, 2.0, 0.25)


def test_WS_03_一维估算():
    结果, 流体, (T烟, T水) = _估算()
    C烟 = 流体['密度'] * 2.0 * 0.25 * 流体['比热']
    C水 = 0.05 * 4400
    assert C烟 * (650 - T烟) == pytest.approx(C水 * (T水 - 450))
    assert 450 < T水 < 650 and 450 < T烟 < 650
    _, _, (T烟顺, T水顺) = _估算({'水流动': '顺流'})
    assert T烟顺 > T烟 and T水顺 < T水               # 逆流换热多
    # Cr = 1 的逆流: ε = NTU/(1+NTU)
    水流量 = C烟 / 4400
    _, _, (T烟1, _) = _估算({'水质量流量': 水流量})
    NTU = 体积换热系数(结果['多孔']['换热']) * 0.5 / C烟
    assert (650 - T烟1) / 200 == pytest.approx(NTU / (1 + NTU))
    # NTU 很大的逆流: 热容流量小的一侧（这里是水，220 < C烟）温差用尽，水出口逼近烟气入口温度
    assert C水 < C烟
    _, _, (_, T水大) = _估算({'管束': dict(管束, 单位体积换热面积=1e5)})
    assert T水大 == pytest.approx(650, abs=0.5)


@pytest.fixture(scope='module')
def 管束床(tmp_path_factory):
    目录 = tmp_path_factory.mktemp('tubebank')
    step = 写STEP(目录 / 'bank.step', lambda occ: [occ.addBox(x0, 0, 0, L, 0.5, 0.5)
                                                     for x0, L in ((0, 1), (1, 2), (3, 1))])
    return step, 目录


def test_WS_04_写算例(管束床, tmp_path):
    step, _ = 管束床
    边界 = generate_boundary(step, 工况(), tmp_path / 'b')
    validate_boundary(边界)
    网格 = generate_mesh(step, 边界, tmp_path / 'm', 0.15)
    信息 = generate_case(step, 边界, 网格, tmp_path / 'c', 步数=20)
    水侧 = 信息['水侧']
    h = 管束换热系数(管束)['h_av']
    assert 水侧['h_av'] == pytest.approx(h) and 水侧['水热容流量'] == pytest.approx(220)
    assert 水侧['水流动'] == '逆流' and 水侧['分段数'] == 40
    assert set(水侧['热阻']) == {'烟气侧', '污垢', '管壁', '水侧'}
    结果 = 校验工况(工况(), True)
    流体 = 流体参数(结果)
    T烟, T水 = 水侧估算(结果, 流体, 2.0, 2.0, 0.25)
    assert 信息['出口温度估算'] == pytest.approx(T烟) and 水侧['水出口温度估算'] == pytest.approx(T水)
    assert 信息['时间步'] == 伪时间步(结果, 流体, 信息['停留时间'])
    k = h / (流体['密度'] * 流体['比热'])
    assert 信息['时间步'] <= 1 / k
    控制 = (tmp_path / 'c/system/controlDict').read_text(encoding='utf-8')
    assert 'type            scalarCodedSource;' in 控制 and 'scalarSemiImplicitSource' not in 控制
    for 片段 in ('const scalar hav = %.12g;' % h, 'const scalar Cw = 220;', 'const scalar TwIn = 450;',
                 'const label N = 40;', 'const bool counter = true;', 'const vector n(1, 0, 0);',
                 'const scalar s0 = 1;', 'const scalar L = 2;'):
        assert 片段 in 控制, 片段
    assert '@' not in 控制
    运行脚本 = (tmp_path / 'c/Allrun').read_text(encoding='utf-8')
    assert 'setpriv --reuid=65534' in 运行脚本 and '"${run_as[@]}" simpleFoam > log.simpleFoam' in 运行脚本
    顺 = generate_case(step, generate_boundary(step, 工况({'水流动': '顺流'}), tmp_path / 'b2'), 网格,
                       tmp_path / 'c2', 步数=20)
    assert 'const bool counter = false;' in (tmp_path / 'c2/system/controlDict').read_text(encoding='utf-8')
    assert 顺['水侧']['水流动'] == '顺流'


# ---- 第 7、8 步造数 ----
水侧信息 = {'水入口温度': 400.0, '水质量流量': 0.5, '水比热': 4400.0, '水热容流量': 2200.0, '水流动': '逆流',
            '分段数': 40, 'h_av': 100.0, '水饱和温度': None, '水出口温度估算': 450.0,
            '监测文件': 'postProcessing/waterSide.dat', '水温分布文件': 'postProcessing/waterProfile.dat'}
七项 = 造数.四项 + ['水侧稳定', '热平衡', '水未沸腾']
十项 = 七项 + ['压差与一维解一致', '出口温度与一维解一致', '水出口温度与一维解一致']


def 水侧算例(root, 水=None, 行数=100, 出口=450.0, 水吸热=110000.0, 烟放热=None, 分布最高=450.0, 改=None, **信息):
    case = 造数.make_case(root, info=dict({'水侧': dict(水侧信息, **(水 or {}))}, **信息))
    pp = case / 'postProcessing'
    行 = ['%d %r %r %r' % (i, 出口, 水吸热, 水吸热 if 烟放热 is None else 烟放热) for i in range(1, 行数 + 1)]
    for i, 值 in (改 or {}).items():
        行[i] = '%d %r %r %r' % (i, 值, 水吸热, 水吸热)
    (pp / 'waterSide.dat').write_text('# Time TwOut Qwater Qgas\n' + '\n'.join(行) + '\n', encoding='utf-8')
    分布 = ['%r %r' % (1 + j / 20, 450 - 50 * j / 20 + (分布最高 - 450) * (j == 0)) for j in range(21)]
    (pp / 'waterProfile.dat').write_text('# s Tw\n' + '\n'.join(分布) + '\n', encoding='utf-8')
    return case


def test_WS_05_检查结果(tmp_path):
    结果 = check_result(水侧算例(tmp_path / 'a'))
    assert [项['项目'] for 项 in 结果['检查']] == 七项 and 结果['通过'] is True, report(结果)
    assert 结果['水出口温度'] == 450.0 and 结果['水吸热量'] == 110000.0 and 结果['水温最高'] == 450.0
    结果 = check_result(水侧算例(tmp_path / 'b'), 对照一维=True)
    assert [项['项目'] for 项 in 结果['检查']] == 十项 and 结果['通过'] is True, report(结果)
    assert json.loads((Path(结果['结果路径'])).read_text(encoding='utf-8'))['水出口温度'] == 450.0

    def 不通过(case, 对照=False):
        return 造数.failed(check_result(case, 对照一维=对照))

    assert 不通过(水侧算例(tmp_path / 'c', 水吸热=100000.0)) == ['热平衡']              # 差 10%
    assert 不通过(水侧算例(tmp_path / 'c2', 水吸热=109000.0)) == []                    # 差 0.9%
    assert 不通过(水侧算例(tmp_path / 'd', 改={-51: 450.02})) == ['水侧稳定']
    assert 不通过(水侧算例(tmp_path / 'e', 烟放热=110500.0)) == ['水侧稳定']           # 水温还在变
    assert 不通过(水侧算例(tmp_path / 'f', 出口=452.0), True) == ['水出口温度与一维解一致']  # 允许 1.5 K
    assert 不通过(水侧算例(tmp_path / 'g', 水={'水饱和温度': 455.0}, 分布最高=456.0)) == ['水未沸腾']
    assert 不通过(水侧算例(tmp_path / 'h', 水={'水饱和温度': 460.0}, 分布最高=456.0)) == []
    非绝热 = 水侧算例(tmp_path / 'i', 水吸热=50000.0, 壁面换热={'方式': '定壁温', '壁温': 300.0})
    结果 = check_result(非绝热)
    assert 造数.failed(结果) == [] and '不对照' in 造数.note(结果, '热平衡')
    缺 = 水侧算例(tmp_path / 'j')
    (缺 / 'postProcessing/waterSide.dat').unlink()
    with pytest.raises(Exception, match='waterSide.dat'):
        check_result(缺)


def test_WS_06_判定(tmp_path):
    for 名, 对照 in (('a', False), ('b', True)):
        case = 水侧算例(tmp_path / 名)
        check_result(case, 对照一维=对照)
        判定 = judge_result(case)
        assert 判定['结论'] == '通过' and 判定['对照一维'] is 对照
    case = 水侧算例(tmp_path / 'c', 水={'水饱和温度': 455.0}, 分布最高=456.0)
    check_result(case)
    判定 = judge_result(case)
    assert 判定['结论'] == '返工' and 判定['回到步骤'] == [2]
    结果 = json.loads((case / '结果.json').read_text(encoding='utf-8'))
    结果['检查'] = 结果['检查'][:5]
    (case / '结果.json').write_text(json.dumps(结果, ensure_ascii=False), encoding='utf-8')
    with pytest.raises(判定错误, match='水侧耦合'):
        judge_result(case)


def _实算(管束床, 名, 工况值):
    from run_case import run_case
    step, 目录 = 管束床
    边界 = generate_boundary(step, 工况值, 目录 / (名 + '-b'))
    validate_boundary(边界)
    网格 = generate_mesh(step, 边界, 目录 / (名 + '-m'), 0.05)
    信息 = generate_case(step, 边界, 网格, 目录 / (名 + '-c'))
    return 信息, check_result(run_case(信息), 对照一维=True)


def test_WS_07_管束逆流实算(能求解, 管束床):
    信息, 结果 = _实算(管束床, '逆流', 工况())
    assert 结果['通过'] is True, report(结果)
    assert [项['项目'] for 项 in 结果['检查']] == 十项
    assert 450 < 结果['水出口温度'] < 650 and 结果['吸热量'] > 0
    assert 结果['水吸热量'] == pytest.approx(结果['吸热量'], rel=0.01)
    test_WS_07_管束逆流实算.水出口 = 结果['水出口温度']


def test_WS_08_顺流直接给总系数(能求解, 管束床):
    值 = 工况({'水流动': '顺流'})
    换 = 值['多孔']['换热']
    换['h_av'] = 管束换热系数(换.pop('管束'))['h_av']
    信息, 结果 = _实算(管束床, '顺流', 值)
    assert 结果['通过'] is True, report(结果)
    assert 结果['水出口温度'] < getattr(test_WS_07_管束逆流实算, '水出口', math.inf)
