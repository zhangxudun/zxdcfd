'''docs/07-检查结果.md 的“OF”类测试点: 真的从第 2 步跑到第 6 步，再交给第 7 步。

T7-09 至 T7-13、T7-23、T7-33、T7-35。几何 A、C、E、G 与 docs/05-写算例.md 2.1 相同，工况是 2.2 的标准工况，
只有 A、G 的壁面改成无滑移，C 的网格尺寸取 0.05，原因见下面 几何、算例表 的注释和 docs/07 2.2。
步数用第 5 步的默认值（T7-23 除外），求解在容器里跑，比较慢。每个算例只求解一次，各测试共用。
'''
import json
import os
import shutil
import subprocess
from pathlib import Path

import gmsh
import pytest

from check_result import check_result
from run_case import run_case

W1 = {'方式': '定热源', '体积热源': -50000.0}
W2 = {'方式': '温差', 'h_av': 2000.0, '壁温顶部': 500.0, '壁温底部': 450.0}
_轴 = {'x': 0, 'y': 1, 'z': 2}


def _长方体(occ):
    occ.addBox(0, 0, 0, 1, 2, 3)


def _叠放(occ):
    for z in range(3):
        occ.addBox(0, 0, z, 1, 2, 1)


def _弯头(occ):
    角 = [(1, 0), (4, 0), (4, 1), (2, 1), (2, 3), (1, 3)]
    点 = [occ.addPoint(x, 0, z) for x, z in 角]
    边 = [occ.addLine(点[i], 点[(i + 1) % len(点)]) for i in range(len(点))]
    截面 = occ.addPlaneSurface([occ.addCurveLoop(边)])
    occ.extrude([(2, 截面)], 0, 1, 0)


def _横放圆柱(occ):
    occ.addCylinder(1, 0, 0.5, 0, 2, 0, 0.5)


# 几何名: (构造, 标注 {(轴, 面中心坐标): 类型}, 网格尺寸)
几何 = {'A': (_长方体, None, None),
        # C 要对照一维压差。默认尺寸（约 0.125）下压差比估算低 3.8%，超过 3% 的门槛: 多孔区边界上
        # 阻力突变，一阶离散把压力梯度抹到了相邻单元，单元中心的速度偏低，压降偏小。偏差随网格加密
        # 单调减小（0.1: 3.4%，0.07: 2.7%，0.05: 2.0%），取 0.05 留出余量。见 docs/07 第 6 节。
        'C': (_叠放, None, 0.05),
        'E': (_弯头, {('z', 0): 'walls', ('x', 4): 'outlet'}, None),
        'G': (_横放圆柱, {('y', 0): 'inlet', ('y', 2): 'outlet'}, 0.05)}
# 算例名: (几何名, 换热, 壁面)
# A、G 是直通道，壁面滑移时解就是均匀流，几十步后速度、压力已经精确，归一化残差只剩舍入噪声（0.1 量级），
# 每步 p 都解满 1000 次，判不了“残差收敛”，也很慢。改成无滑移，流场有边界层，残差能正常下降。
算例表 = {'A': ('A', None, '无滑移'), 'C1': ('C', W1, '滑移'), 'C2': ('C', W2, '滑移'),
          'E': ('E', None, '滑移'), 'G': ('G', None, '无滑移')}


def 工况(换热=None, 壁面='滑移'):
    值 = {'入口': {'速度': 5.0, '温度': 600.0},
          '出口': {'表压': -200.0},
          '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                   '湍流粘度': 0.002, '湍流普朗特数': 0.85},
          '壁面': 壁面}
    if 换热:
        值['多孔'] = {'粘性阻力系数': 1000.0, '惯性阻力系数': 20.0, '横向倍数': 10.0, '换热': 换热}
    return 值


def 写算例(目录, 名, 步数=None):
    '''第 2 至 5 步，返回 (第 5 步的返回值, 第 4 步的网格统计)。'''
    from generate_boundary import generate_boundary, 识别几何
    from generate_case import generate_case
    from generate_mesh import generate_mesh
    from generate_structure import _写出_单位m
    from validate_boundary import validate_boundary
    几何名, 换热, 壁面 = 算例表[名]
    构造, 改, 网格尺寸 = 几何[几何名]
    模型 = os.path.join(str(目录), 几何名 + '.step')
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add('geo')
        构造(gmsh.model.occ)
        gmsh.model.occ.synchronize()
        _写出_单位m(模型)
    finally:
        gmsh.finalize()
    标注 = None
    if 改:
        面 = 识别几何(模型)['面']
        标注 = {'面': {}}
        for (轴, 值), 类 in 改.items():
            命中 = [f['编号'] for f in 面 if abs(f['中心'][_轴[轴]] - 值) <= 1e-3]
            assert len(命中) == 1, '中心在 %s=%s 的面应恰好一个，实际 %d 个' % (轴, 值, len(命中))
            标注['面'][命中[0]] = 类
    边界 = generate_boundary(模型, 工况(换热, 壁面), str(目录), 标注)
    validate_boundary(边界)
    网格 = generate_mesh(模型, 边界, str(目录), 网格尺寸)
    参数 = {} if 步数 is None else {'步数': 步数}
    return generate_case(模型, 边界, 网格, Path(str(目录)) / '算例', **参数), 网格


def docker_ok():
    if shutil.which('docker') is None:
        return False
    return subprocess.run(['docker', 'info'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


@pytest.fixture(scope='module')
def solved(tmp_path_factory):
    '''solved(算例名, 步数=None) 返回 {'信息': 第 5 步返回值, '网格': 第 4 步统计, '求解': 第 6 步返回值}。'''
    if not docker_ok():
        if os.environ.get('CI'):
            pytest.fail('CI 里不允许跳过 OF 类测试: 不能运行 docker')
        pytest.skip('不能运行 docker，OF 类测试跳过（跳过不算通过）')
    done = {}

    def get(name, steps=None):
        key = (name, steps)
        if key not in done:
            try:
                info, mesh = 写算例(tmp_path_factory.mktemp('%s-%s' % key), name, steps)
                done[key] = {'信息': info, '网格': mesh, '求解': run_case(info)}
            except Exception as error:          # 失败也记下来，同一个算例不重复求解
                done[key] = error
        if isinstance(done[key], Exception):
            raise done[key]
        return done[key]
    return get


def report(result):
    return json.dumps(result, ensure_ascii=False, indent=1)


def case_files(case):
    return {str(p.relative_to(case)): p.read_bytes() for p in Path(case).rglob('*')
            if p.is_file() and p.name != '结果.json'}


def last_row(path):
    rows = [line.split() for line in Path(path).read_text(encoding='utf-8').splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    return float(rows[-1][-1]), len(rows)


@pytest.mark.parametrize('name', ['C1', 'C2', 'A', 'E', 'G'])
def test_real_case_passes(solved, name):                                        # T7-09 至 T7-12
    data = solved(name)
    compare = name in ('C1', 'C2')
    result = check_result(data['求解'], 对照一维=compare)
    assert len(result['检查']) == (6 if compare else 4)
    assert result['通过'] is True, report(result)
    if name == 'C1':
        volume = sum(data['网格']['体'][t]['网格体积'] for t in data['网格']['多孔体'])
        assert result['吸热量'] == pytest.approx(50000 * volume, rel=0.03), report(result)
    if name == 'A':
        assert abs(result['入口流量'] - 10) < 0.1, report(result)


def test_numbers_come_from_monitor_files(solved):                               # T7-13
    run = solved('C1')['求解']
    result = check_result(run)
    value = {name: last_row(path) for name, path in run['监测数据'].items()}
    assert result['迭代步数'] == run['完成步数'] == value['inletP'][1]
    assert result['入口流量'] == -value['inletFlow'][0]
    assert result['出口流量'] == value['outletFlow'][0]
    assert result['出口温度'] == value['outletT'][0]


def test_return_value_and_path_give_same_file(solved):                          # T7-33
    run = solved('C1')['求解']
    case = Path(run['算例目录'])
    first = check_result(run, 对照一维=True)
    content = (case / '结果.json').read_bytes()
    assert check_result(str(case), 对照一维=True) == first
    assert (case / '结果.json').read_bytes() == content


def test_solver_output_is_untouched(solved):                                    # T7-35
    run = solved('C1')['求解']
    before = case_files(run['算例目录'])
    assert 'log.simpleFoam' in before and 'constant/polyMesh/boundary' in before and '算例信息.json' in before
    check_result(run, 对照一维=True)
    assert case_files(run['算例目录']) == before


def test_twenty_steps_is_not_steady(solved):                                    # T7-23
    result = check_result(solved('C1', 20)['求解'])
    items = {item['项目']: item for item in result['检查']}
    assert list(items) == ['流量守恒', '入口流量', '监测值稳定', '残差收敛'], report(result)
    assert items['监测值稳定']['通过'] is False and '共 20 步' in items['监测值稳定']['说明'], report(result)
    assert all(isinstance(items[name]['通过'], bool) and items[name]['说明'] for name in ('流量守恒', '入口流量', '残差收敛'))
    assert result['通过'] is False and result['迭代步数'] == 20
