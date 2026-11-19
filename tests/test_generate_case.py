import json
import re
from pathlib import Path

import pytest

from generate_boundary import generate_boundary
from generate_case import 算例错误, generate_case
from generate_mesh import generate_mesh
from generate_structure import generate_step, 默认参数
from validate_boundary import validate_boundary


@pytest.fixture(scope='module', params=['S1', 'S2'])
def case_input(request, tmp_path_factory):
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
    return request.param, step, boundary, mesh


def snapshot(root):
    return {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_generate_case_from_real_mesh(case_input, tmp_path):
    shape, step, boundary, mesh = case_input
    out = tmp_path / 'case'
    result = generate_case(step, boundary, mesh, out, 步数=200)
    files = snapshot(out)
    expected = {'inlet', 'outlet', 'walls'} | ({'symmetry'} if shape == 'S2' else set())
    assert set(result['边界']) == set(mesh['边界']) == expected
    for name in ('U', 'p', 'T'):
        text = (out / '0' / name).read_text(encoding='utf-8')
        assert set(re.findall(r'^    (\w+)\s*\n    \{', text, re.M)) == expected
    assert 'refValue uniform -5;' in (out / '0/U').read_text(encoding='utf-8')
    control = (out / 'system/controlDict').read_text(encoding='utf-8')
    assert re.search(r'endTime\s+200;', control)
    assert 'scalarTransport' in control and 'porous_zone' in control
    assert 'porous_zone' in (out / 'system/fvOptions').read_text(encoding='utf-8')
    assert (out / 'mesh.msh').read_bytes() == Path(mesh['网格路径']).read_bytes()
    summary = json.loads((out / '算例信息.json').read_text(encoding='utf-8'))
    assert summary == {k: v for k, v in result.items() if k not in ('算例目录', '算例信息路径')}
    assert summary['多孔区'] == 'porous_zone'
    assert summary['多孔区流速'] == pytest.approx(5 if shape == 'S1' else 1.25)
    assert summary['多孔区长度'] == pytest.approx(3.2)
    assert summary['出口温度估算'] == pytest.approx(550.6112 if shape == 'S1' else 402.4450, abs=0.001)
    generate_case(step, boundary['标注路径'], mesh['网格路径'], out, 工况=boundary['工况'], 步数=200)
    assert snapshot(out) == files


@pytest.mark.parametrize('steps', [0, -1, 2.5, '200', True, float('nan')])
def test_invalid_steps_preserve_existing_case(case_input, tmp_path, steps):
    _, step, boundary, mesh = case_input
    out = tmp_path / 'case'
    out.mkdir()
    (out / 'Allrun').write_text('old', encoding='utf-8')
    old = snapshot(out)
    with pytest.raises(算例错误, match='步数'):
        generate_case(step, boundary, mesh, out, 步数=steps)
    assert snapshot(out) == old
