import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest

from run_case import run_case, 报错关键字, 求解错误

监测量 = ['inletFlow', 'inletP', 'outletFlow', 'outletT']
边界 = {'inlet': 'patch', 'outlet': 'patch', 'symmetry': 'symmetry', 'walls': 'wall'}
边界文件 = '''FoamFile
{
    version     2.0;
    class       polyBoundaryMesh;
    object      boundary;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

4
(
    walls
    {
        type wall;
        physicalType    patch;
        nFaces          2142;
    }
    symmetry
    {
        type symmetry;
        physicalType    patch;
        nFaces          732;
    }
    inlet
    {
        type            patch;
        physicalType    patch;
        nFaces          152;
    }
    outlet
    {
        type            patch;
        physicalType    patch;
        nFaces          88;
    }
)
'''


def stub(steps=20, converged=False, extra=''):
    '''假的 Allrun: 写出三个合格的日志和四个监测量的数据，extra 是追加的 bash 命令。'''
    return '''#!/bin/bash
cd "$(dirname "$0")"
printf 'tet  :10\\n--> FOAM Warning : \\n    Found 4 undefined faces in mesh; adding to default patch defaultFaces\\nEnd\\n' > log.gmshToFoam
printf 'Mesh OK.\\n\\nEnd\\n' > log.checkMesh
{ for i in $(seq 1 %d); do echo "Time = $i"; echo; done; %s echo End; } > log.simpleFoam
for m in %s; do
    mkdir -p postProcessing/$m/0
    { echo '# Time value'; for i in $(seq 1 %d); do echo "$i 1.5"; done; } > postProcessing/$m/0/surfaceFieldValue.dat
done
%s
''' % (steps, 'echo "SIMPLE solution converged in %d iterations";' % steps if converged else '',
       ' '.join(监测量), steps, extra)


def make_case(root, script, name='case'):
    case = root / name
    (case / '0').mkdir(parents=True)
    (case / '0' / 'U').write_text('U', encoding='utf-8')
    (case / 'Allrun').write_text(script, encoding='utf-8')
    (case / 'mesh.msh').write_text('$MeshFormat', encoding='utf-8')
    (case / '算例信息.json').write_text(
        json.dumps({'边界': 边界, '监测量': 监测量, '步数': 20}, ensure_ascii=False), encoding='utf-8')
    return case


def snapshot(root, names=('0/U', 'Allrun', 'mesh.msh', '算例信息.json')):
    return {name: (root / name).read_bytes() for name in names}


def logs(case):
    return sorted(p.name for p in case.glob('log.*'))


# ---------- 桩: 不需要 OpenFOAM ----------

@pytest.mark.parametrize('converged', [False, True])
def test_stub_pass(tmp_path, monkeypatch, converged):
    script = stub(12 if converged else 20, converged,
                  "mkdir -p constant/polyMesh; cat > constant/polyMesh/boundary <<'EOF'\n" + 边界文件 + 'EOF')
    case = make_case(tmp_path, script, '算例')
    old = snapshot(case)
    monkeypatch.chdir(tmp_path)
    before = sorted(os.listdir(tmp_path))
    result = run_case('算例', 运行方式='本机')
    assert result['算例目录'] == str(case)
    assert result['完成步数'] == (12 if converged else 20) and result['步数'] == 20
    assert result['提前收敛'] is converged
    assert result['边界'] == 边界
    assert sorted(result['监测数据']) == 监测量
    assert all(os.path.isfile(p) for p in list(result['日志'].values()) + list(result['监测数据'].values()))
    assert logs(case) == ['log.checkMesh', 'log.gmshToFoam', 'log.simpleFoam']
    assert snapshot(case) == old and sorted(os.listdir(tmp_path)) == before
    again = run_case({'算例目录': str(case)}, 运行方式='本机')
    assert again == result


def test_old_results_are_cleared(tmp_path):
    case = make_case(tmp_path, stub())
    (case / 'log.simpleFoam').write_text('--> FOAM FATAL ERROR\n', encoding='utf-8')
    (case / 'log.old').write_text('old', encoding='utf-8')
    for name in ('100', 'postProcessing/gone/0', 'constant/polyMesh', 'notes'):
        (case / name).mkdir(parents=True)
    run_case(case, 运行方式='本机')
    assert not (case / '100').exists() and not (case / 'postProcessing/gone').exists()
    assert not (case / 'constant/polyMesh').exists() and not (case / 'log.old').exists()
    assert (case / '0/U').exists() and (case / 'notes').is_dir()


BAD = [
    ("printf 'Failed 2 mesh checks.\\n\\nEnd\\n' > log.checkMesh", 'checkMesh 不通过.*第 4 步'),
    ('rm log.simpleFoam', '没有 simpleFoam 的日志'),
    (': > log.simpleFoam', 'simpleFoam 的日志是空的'),
    ('echo tail >> log.simpleFoam', 'simpleFoam 的日志.*不以 End 结束'),
    ('rm -r postProcessing/outletT', '没有监测量 outletT'),
    ("echo '# Time value' > postProcessing/inletP/0/surfaceFieldValue.dat", '监测量 inletP 没有数据行'),
    ("sed -i '$ s/1.5/nan/' postProcessing/outletFlow/0/surfaceFieldValue.dat", '监测量 outletFlow 第 20 行'),
    ("sed -i '$ d' postProcessing/inletFlow/0/surfaceFieldValue.dat", '监测量 inletFlow 有 19 行'),
    ("mkdir -p constant/polyMesh; printf 'inlet\\n{\\n type patch;\\n}\\ndefaultFaces\\n{\\n type patch;\\n}\\n'"
     ' > constant/polyMesh/boundary', '边界与算例信息不一致.*defaultFaces'),
] + [("{ echo '--> %s'; cat log.simpleFoam; } > t; mv t log.simpleFoam" % word,
      'OpenFOAM' if word == 'command not found' else 'simpleFoam 的日志.*' + word) for word in 报错关键字]


@pytest.mark.parametrize('extra, match', BAD)
def test_stub_exit_zero_but_not_ok(tmp_path, extra, match):
    case = make_case(tmp_path, stub(extra=extra))
    old = snapshot(case)
    with pytest.raises(求解错误, match=match):
        run_case(case, 运行方式='本机')
    assert snapshot(case) == old


@pytest.mark.parametrize('script, match, kept', [
    ('echo broken > log.gmshToFoam\nexit 1\n', '停在 gmshToFoam，退出码 1', ['log.gmshToFoam']),
    ('echo "Allrun: line 9: gmshToFoam: command not found" > log.gmshToFoam\nexit 127\n',
     '停在 gmshToFoam，退出码 127。环境里找不到 OpenFOAM', ['log.gmshToFoam']),
    ('echo End > log.gmshToFoam\necho "Mesh OK." > log.checkMesh\n'
     'printf -- "--> FOAM FATAL IO ERROR\\nkeyword inlet is undefined\\n" > log.simpleFoam\nexit 1\n',
     '停在 simpleFoam，退出码 1.*keyword inlet is undefined',
     ['log.checkMesh', 'log.gmshToFoam', 'log.simpleFoam']),
    ('echo oops\nexit 3\n', 'Allrun 没有跑到 gmshToFoam，退出码 3: oops', []),
])
def test_stub_nonzero_exit(tmp_path, script, match, kept):
    case = make_case(tmp_path, '#!/bin/bash\ncd "$(dirname "$0")"\n' + script)
    with pytest.raises(求解错误, match=match):
        run_case(case, 运行方式='本机')
    assert logs(case) == kept


def test_timeout_kills_the_run(tmp_path):
    case = make_case(tmp_path, '#!/bin/bash\ncd "$(dirname "$0")"\necho $$ > pid\nsleep 60 &\necho $! >> pid\nwait\n')
    start = time.time()
    with pytest.raises(求解错误, match='超时'):
        run_case(case, 时限=2, 运行方式='本机')
    assert time.time() - start < 20
    time.sleep(0.5)
    for pid in (case / 'pid').read_text().split():
        assert not os.path.exists('/proc/%s' % pid) or 'Z' in open('/proc/%s/stat' % pid).read().split(')')[-1][:3]


@pytest.mark.parametrize('missing', ['目录', 'Allrun', 'mesh.msh', '算例信息.json'])
def test_incomplete_case_is_not_run(tmp_path, missing):
    case = make_case(tmp_path, stub())
    if missing == '目录':
        shutil.rmtree(case)
    else:
        (case / missing).unlink()
    with pytest.raises(求解错误, match='找不到算例目录' if missing == '目录' else missing):
        run_case(case, 运行方式='本机')
    assert not case.exists() or logs(case) == []


@pytest.mark.parametrize('kwargs, match', [
    ({'时限': 0}, '时限'), ({'时限': True}, '时限'), ({'时限': float('nan')}, '时限'),
    ({'时限': '60'}, '时限'), ({'运行方式': '远程'}, '运行方式'),
])
def test_invalid_arguments(tmp_path, kwargs, match):
    case = make_case(tmp_path, stub())
    with pytest.raises(求解错误, match=match):
        run_case(case, **dict({'运行方式': '本机'}, **kwargs))
    assert logs(case) == []


def test_container_mode_without_docker(tmp_path, monkeypatch):
    case = make_case(tmp_path, stub())
    (case / 'log.old').write_text('old', encoding='utf-8')
    monkeypatch.setenv('PATH', str(tmp_path))
    with pytest.raises(求解错误, match='找不到 docker'):
        run_case(case)
    assert logs(case) == ['log.old']


# ---------- OF: 用容器真的运行 OpenFOAM ----------

def docker_ok():
    if shutil.which('docker') is None:
        return False
    return subprocess.run(['docker', 'info'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


@pytest.fixture(scope='module', params=['S1', 'S2'])
def real_case(request, tmp_path_factory):
    if not docker_ok():
        if os.environ.get('CI'):
            pytest.fail('CI 里不允许跳过 OF 类测试: 不能运行 docker')
        pytest.skip('不能运行 docker，OF 类测试跳过（跳过不算通过）')
    from generate_boundary import generate_boundary
    from generate_case import generate_case
    from generate_mesh import generate_mesh
    from generate_structure import generate_step, 默认参数
    from validate_boundary import validate_boundary
    root = tmp_path_factory.mktemp(request.param)
    step = generate_step(request.param, 默认参数(request.param), root, 'duct')
    condition = {
        '入口': {'速度': 5.0, '温度': 600.0},
        '出口': {'表压': -200.0},
        '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                 '湍流粘度': 0.002, '湍流普朗特数': 0.85},
        '壁面': '滑移',
        '多孔': {'粘性阻力系数': 1000.0, '惯性阻力系数': 20.0, '横向倍数': 10.0,
                 '换热': {'方式': '定热源', '体积热源': -50000.0}},
    }
    boundary = generate_boundary(step, condition, root)
    validate_boundary(boundary)
    mesh = generate_mesh(step, boundary, root)
    return request.param, generate_case(step, boundary, mesh, root / '算例', 步数=20), mesh


def case_files(case):
    return {str(p.relative_to(case)): p.read_bytes() for p in case.rglob('*') if p.is_file()}


def test_real_case_runs(real_case):
    shape, info, mesh = real_case
    case = Path(info['算例目录'])
    old = case_files(case)
    result = run_case(info)
    expected = {'inlet': 'patch', 'outlet': 'patch', 'walls': 'wall'}
    if shape == 'S2':
        expected['symmetry'] = 'symmetry'
    assert result['边界'] == info['边界'] == expected
    assert result['完成步数'] == 20 and result['提前收敛'] is False
    assert sorted(result['监测数据']) == 监测量
    assert all(os.path.isfile(p) for p in result['日志'].values())
    gmsh_log = (case / 'log.gmshToFoam').read_text()
    assert 'tet  :%d' % mesh['单元数'] in gmsh_log
    assert 'porous_zone' in (case / 'constant/polyMesh/cellZones').read_text()
    assert all((case / '20' / name).is_file() for name in ('U', 'p', 'T', 'phi'))
    new = case_files(case)
    assert {name: new[name] for name in old} == old
    added = {name.split('/')[0] if not name.startswith('constant/') else 'constant/polyMesh'
             for name in set(new) - set(old)}
    assert added == {'log.gmshToFoam', 'log.checkMesh', 'log.simpleFoam', 'postProcessing',
                     'constant/polyMesh', '20'}
    assert all(os.access(p, os.R_OK | os.W_OK) for p in case.rglob('*'))
    again = run_case(case)
    assert again == result


def test_real_broken_case_stops_at_simplefoam(real_case, tmp_path):
    _, info, _ = real_case
    case = tmp_path / 'case'
    shutil.copytree(info['算例目录'], case, ignore=shutil.ignore_patterns(
        'log.*', 'postProcessing', 'polyMesh', '[1-9]*'))
    (case / 'system/fvSchemes').unlink()
    with pytest.raises(求解错误, match='停在 simpleFoam'):
        run_case(case)
    assert (case / 'log.gmshToFoam').is_file()
