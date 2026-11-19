import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import judge_result as module
from judge_result import _命令行, judge_result, 判定错误

项目 = ['流量守恒', '入口流量', '监测值稳定', '压差与一维解一致', '出口烟温与一维解一致']
步骤 = {'流量守恒': [5, 4], '监测值稳定': [5, 4], '入口流量': [2, 4],
        '压差与一维解一致': [2, 4], '出口烟温与一维解一致': [2, 4]}
键 = ['结论', '回到步骤', '对照一维', '原因']


def r0(fail=(), n=5):
    '''docs/08 的基准结果 R0；fail 里的检查项不通过，n=3 是“R0 三项”。'''
    return {'通过': not fail, '迭代步数': 100, '进出口压差': 150.0, '入口流量': 10.0, '出口流量': 10.0,
            '出口温度': 580.0, '吸热量': 110000.0, '多孔压降估算': 150.0, '出口温度估算': 580.0,
            '检查': [{'项目': 名, '通过': 名 not in fail, '说明': 名 + '的说明'} for 名 in 项目[:n]]}


def write(root, result, name='case'):
    case = root / name
    case.mkdir(parents=True, exist_ok=True)
    text = result if isinstance(result, str) else json.dumps(result, ensure_ascii=False)
    (case / '结果.json').write_text(text, encoding='utf-8')
    return case


def files(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob('*') if p.is_file()}


def raised(func, *args, **kwargs):
    try:
        func(*args, **kwargs)
    except Exception as e:
        return e
    pytest.fail('%s 没有报错' % func.__name__)


# ---------- 造数: 通过 ----------

@pytest.mark.parametrize('n', [5, 3])
def test_pass(tmp_path, n):                                   # T8-01、T8-02、T8-03
    result = dict(r0(n=n), 算例目录=str(tmp_path), 结果路径=str(tmp_path / '结果.json'))
    verdict = judge_result(result)
    assert verdict == {'结论': '通过', '回到步骤': [], '对照一维': n == 5, '原因': []}
    assert judge_result(write(tmp_path, r0(n=n))) == verdict
    assert judge_result(str(tmp_path / 'case')) == verdict


# ---------- 造数: 检查不通过 ----------

@pytest.mark.parametrize('name, n, words', [
    ('流量守恒', 5, ['步数', '加密网格']), ('监测值稳定', 5, ['步数', '加密网格']),
    ('入口流量', 5, ['入口面', '标注']), ('入口流量', 3, ['入口面', '标注']),
    ('压差与一维解一致', 5, ['多孔体', '流向', '一维估算是否适用']),
    ('出口烟温与一维解一致', 5, ['多孔体', '流向', '一维估算是否适用']),
])
def test_one_failed_check(name, n, words):                    # T8-06 至 T8-09、T8-14
    verdict = judge_result(r0([name], n))
    assert verdict['结论'] == '返工' and verdict['对照一维'] is (n == 5)
    assert verdict['回到步骤'] == sorted(步骤[name])
    assert len(verdict['原因']) == 1
    reason = verdict['原因'][0]
    assert name in reason['现象'] and reason['回到步骤'] == 步骤[name]
    assert all(word in reason['处理'] for word in words)
    assert reason['说明'] == name + '的说明'


@pytest.mark.parametrize('fail, steps', [
    (['流量守恒', '监测值稳定'], [4, 5]),
    (['流量守恒', '压差与一维解一致'], [2, 4, 5]),
    (项目, [2, 4, 5]),
])
def test_several_failed_checks(fail, steps):                  # T8-10、T8-11、T8-12
    verdict = judge_result(r0(fail))
    assert verdict['结论'] == '返工' and verdict['回到步骤'] == steps
    assert [next(名 for 名 in 项目 if 名 in r['现象']) for r in verdict['原因']] == fail
    assert [r['回到步骤'] for r in verdict['原因']] == [步骤[名] for 名 in fail]


def test_note_is_copied():                                    # T8-13
    result = r0(['流量守恒'])
    result['检查'][0]['说明'] = '入口流量 10 m³/s，出口流量 10.011 m³/s，相差 0.011（允许 0.01）'
    assert judge_result(result)['原因'][0]['说明'] == result['检查'][0]['说明']


# ---------- 触发: 前面的步骤报错 ----------

好日志 = ("echo End > log.gmshToFoam\nprintf 'Mesh OK.\\n' > log.checkMesh\n"
          "printf 'Time = 1\\nEnd\\n' > log.simpleFoam\n")


def of_case(root, script):
    '''只够 run_case 走到判定的算例目录；script 是 Allrun 的内容，None 表示没有 Allrun。'''
    case = root / 'of'
    case.mkdir()
    if script is not None:
        (case / 'Allrun').write_text('#!/bin/bash\ncd "$(dirname "$0")"\n' + script, encoding='utf-8')
    (case / 'mesh.msh').write_text('$MeshFormat', encoding='utf-8')
    (case / '算例信息.json').write_text(
        json.dumps({'边界': {}, '监测量': [], '步数': 1}, ensure_ascii=False), encoding='utf-8')
    return case


def d0(root, rows=100, post=True, **change):
    '''docs/07 的基准数据 D0；change 里值为 None 的键改成空，值为 ... 的键去掉。'''
    case = root / 'd0'
    info = {'密度': 0.5, '比热': 1100.0, '入口温度': 600.0, '出口表压': -200.0, '理论入口体积流量': 10.0,
            '多孔压降估算': 150.0, '出口温度估算': 580.0, '步数': rows,
            '监测量': ['inletFlow', 'inletP', 'outletFlow', 'outletT']}
    info.update(change)
    info = {k: v for k, v in info.items() if v is not ...}
    case.mkdir()
    (case / '算例信息.json').write_text(json.dumps(info, ensure_ascii=False), encoding='utf-8')
    if post:
        for name, value in {'inletP': -100, 'inletFlow': -10, 'outletFlow': 10, 'outletT': 580}.items():
            folder = case / 'postProcessing' / name / '0'
            folder.mkdir(parents=True)
            lines = ['# Region', '# Time value'] + ['%d %s' % (i, value) for i in range(1, rows + 1)]
            (folder / 'surfaceFieldValue.dat').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return case


def f1(root, monkeypatch):
    from validate_boundary import validate_boundary
    return raised(validate_boundary, {'面': {'1': 'walls'}})


def f2(root, monkeypatch):
    from generate_mesh import generate_mesh
    return raised(generate_mesh, root / 'x.step', {'面': {'1': '不支持'}}, root / 'out')


def f3(root, monkeypatch):
    from generate_case import generate_case
    return raised(generate_case, root / 'x.step', {'面': {'1': 'inlet'}}, root / 'mesh.msh', root / 'out')


def f4(root, monkeypatch):
    from generate_mesh import generate_mesh
    return raised(generate_mesh, root / 'x.step', {'面': {'1': 'inlet'}}, root / 'out', -1)


def f5(root, monkeypatch):
    from generate_case import generate_case
    return raised(generate_case, root / 'x.step', {'面': {'1': 'inlet'}}, root / 'mesh.msh', root / 'out',
                  步数=0)


def run(root, script, **kwargs):
    from run_case import run_case
    return raised(run_case, of_case(root, script), **dict({'运行方式': '本机'}, **kwargs))


def f6(root, monkeypatch):
    return run(root, 好日志 + "printf 'Failed 2 mesh checks.\\n' > log.checkMesh\n")


def f7_no_openfoam(root, monkeypatch):
    return run(root, "echo 'Allrun: line 9: gmshToFoam: command not found' > log.gmshToFoam\nexit 127\n")


def f7_timeout(root, monkeypatch):
    return run(root, 'sleep 60\n', 时限=1)


def f7_no_docker(root, monkeypatch):
    monkeypatch.setenv('PATH', str(root))
    return run(root, 好日志, 运行方式='容器')


def f8(root, monkeypatch):
    return run(root, None)


def f9(root, monkeypatch):
    return run(root, 好日志 + "printf -- '--> FOAM FATAL ERROR\\nboom\\nEnd\\n' > log.simpleFoam\n")


def f10(root, monkeypatch):
    from check_result import check_result
    return raised(check_result, d0(root, 多孔压降估算=None, 出口温度估算=None), 对照一维=True)


def f11(root, monkeypatch):
    from check_result import check_result
    return raised(check_result, d0(root, 密度=...))


def f11_no_monitor_list(root, monkeypatch):
    from check_result import check_result
    return raised(check_result, d0(root, 监测量=...))


def f12(root, monkeypatch):
    from check_result import check_result
    return raised(check_result, d0(root, post=False))


@pytest.mark.parametrize('trigger, kind, steps, words', [
    (f1, '校验错误', [2], ['在标注里指定面编号']),
    (f2, '标注错误', [2], []), (f3, '工况错误', [2], []),
    (f4, '网格错误', [4], ['减小网格尺寸', '修几何']),
    (f5, '算例错误', [5], []),
    (f6, '求解错误', [4], ['减小网格尺寸', '修几何']),
    (f7_no_openfoam, '求解错误', [6], ['算例不用改']),
    (f7_timeout, '求解错误', [6], ['算例不用改']),
    (f7_no_docker, '求解错误', [6], ['算例不用改']),
    (f8, '求解错误', [5], []), (f9, '求解错误', [5, 4], []),
    (f10, '结果错误', [7], ['不做一维对照', '三项基本检查']),
    (f11, '结果错误', [5], []), (f11_no_monitor_list, '结果错误', [5], []),
    (f12, '结果错误', [6], []),
])
def test_failed_step(tmp_path, monkeypatch, trigger, kind, steps, words):     # T8-16 至 T8-25、T8-32
    error = trigger(tmp_path, monkeypatch)
    assert type(error).__name__ == kind
    args, text, old = error.args, str(error), files(tmp_path)
    verdict = judge_result(error)
    assert verdict['结论'] == '返工' and verdict['对照一维'] is False
    assert verdict['回到步骤'] == sorted(steps) and len(verdict['原因']) == 1
    reason = verdict['原因'][0]
    assert kind in reason['现象'] and reason['回到步骤'] == steps and reason['说明'] == text
    assert all(word in reason['处理'] for word in words)
    assert error.args == args and str(error) == text and files(tmp_path) == old


# ---------- 造数: 第 7 步的真实产物交给第 8 步 ----------

def test_check_result_output_is_accepted(tmp_path):
    from check_result import check_result
    case = d0(tmp_path)
    for compare in (True, False):
        result = check_result(case, 对照一维=compare)
        expected = {'结论': '通过', '回到步骤': [], '对照一维': compare, '原因': []}
        assert judge_result(result) == expected and judge_result(case) == expected


def test_check_result_failure_is_sent_back(tmp_path):
    from check_result import check_result
    result = check_result(d0(tmp_path, rows=20))
    verdict = judge_result(result)
    assert verdict['结论'] == '返工' and verdict['回到步骤'] == [4, 5]
    assert [r['现象'] for r in verdict['原因']] == ['监测值稳定不通过']
    assert verdict['原因'][0]['说明'] == result['检查'][2]['说明']


# ---------- 不能判定 ----------

def test_missing_or_unreadable_result(tmp_path):              # T8-26
    with pytest.raises(判定错误, match='找不到算例目录'):
        judge_result(tmp_path / '没有')
    (tmp_path / 'empty').mkdir()
    with pytest.raises(判定错误, match='没有 结果.json.*第 7 步'):
        judge_result(tmp_path / 'empty')
    for text in ('{', '[1, 2]'):
        with pytest.raises(判定错误, match='结果.json 读不了'):
            judge_result(write(tmp_path, text))


def broken(change):
    result = r0()
    change(result)
    return result


@pytest.mark.parametrize('result, match', [
    (broken(lambda r: r['检查'][1].update(通过=False)), '最外层的 通过 与各项不一致'),                   # T8-27
    (broken(lambda r: r.update(通过=False)), '最外层的 通过 与各项不一致'),
    (broken(lambda r: r.pop('通过')), '没有 通过.*第 7 步'),                                             # T8-28
    (broken(lambda r: r.pop('检查')), '没有 检查.*第 7 步'),
    (broken(lambda r: r.update(检查={'流量守恒': True})), '检查 不是列表.*第 7 步'),
    (broken(lambda r: r.update(通过='真')), '通过 不是布尔值.*第 7 步'),
    (broken(lambda r: r['检查'][0].pop('项目')), '第 1 项没有 项目.*第 7 步'),
    (broken(lambda r: r['检查'][2].pop('通过')), '第 3 项没有 通过.*第 7 步'),
    (broken(lambda r: r['检查'][2].update(通过=1)), '第 3 项.*通过 不是布尔值'),
    (broken(lambda r: r.update(检查=[])), '检查 不是规定的三项或五项.*空'),                              # T8-29
    (broken(lambda r: r.update(检查=r['检查'][:2])), '检查 不是规定的三项或五项.*2 项'),
    (broken(lambda r: r.update(检查=r['检查'][:4])), '检查 不是规定的三项或五项.*4 项'),
    (broken(lambda r: r.update(检查=[r['检查'][i] for i in (1, 0, 2)])), '检查 不是规定的三项或五项.*3 项'),
    (broken(lambda r: r['检查'][4].update(项目='别的检查')), '检查.*不认识的项目 别的检查'),
    (broken(lambda r: r['检查'][1].update(项目='流量守恒')), '检查.*流量守恒 重复'),
])
def test_untrustworthy_result(tmp_path, result, match):
    with pytest.raises(判定错误, match=match):
        judge_result(result)
    with pytest.raises(判定错误, match=match):
        judge_result(write(tmp_path, result))


@pytest.mark.parametrize('value', [None, 5, 1.5, True, [r0()], b'case'])
def test_invalid_input(value):                                # T8-30
    with pytest.raises(判定错误, match='输入必须是'):
        judge_result(value)


@pytest.mark.parametrize('error', [KeyError('通过'), ValueError('普通的错误'), RuntimeError('boom'),
                                   判定错误('不能判定')])
def test_unknown_exception(error):                            # T8-30
    with pytest.raises(判定错误, match='不是第 2–7 步的异常.*' + type(error).__name__):
        judge_result(error)


# ---------- 输入与产物 ----------

def test_nothing_is_changed(tmp_path, monkeypatch):           # T8-31
    cases = [write(tmp_path, r0(), 'ok'), write(tmp_path, r0(['流量守恒']), 'redo'), tmp_path / 'none']
    cases[2].mkdir()
    monkeypatch.chdir(tmp_path)
    old = files(tmp_path)
    assert judge_result(cases[0])['结论'] == '通过' and judge_result(cases[1])['结论'] == '返工'
    with pytest.raises(判定错误):
        judge_result(cases[2])
    for result in (r0(), r0(项目)):
        kept = copy.deepcopy(result)
        judge_result(result)
        assert result == kept
    assert files(tmp_path) == old


def test_verdict_is_json(tmp_path):                           # T8-33
    for result in (r0(), r0(项目)):
        verdict = judge_result(result)
        text = json.dumps(verdict, ensure_ascii=False, allow_nan=False)
        assert json.loads(text) == verdict and '\\u' not in text
        assert list(verdict) == 键
        assert all(sorted(r) == sorted(['现象', '回到步骤', '处理', '说明']) for r in verdict['原因'])
    verdict['原因'][0]['回到步骤'].append(9)
    assert judge_result(r0(项目))['原因'][0]['回到步骤'] == [5, 4]


def test_no_step_is_rerun(tmp_path, monkeypatch):             # T8-34
    import check_result
    import run_case
    import validate_boundary

    def boom(*args, **kwargs):
        raise AssertionError('第 8 步不应该重跑任何步骤')

    monkeypatch.setattr(check_result, 'check_result', boom)
    monkeypatch.setattr(run_case, 'run_case', boom)
    monkeypatch.setattr(validate_boundary, 'validate_boundary', boom)
    for name in ('Popen', 'run', 'call', 'check_call', 'check_output'):
        monkeypatch.setattr(subprocess, name, boom)
    monkeypatch.setattr(os, 'system', boom)
    assert judge_result(write(tmp_path, r0()))['结论'] == '通过'
    assert judge_result(r0(['流量守恒']))['结论'] == '返工'
    steps = {'generate_structure', 'generate_boundary', 'validate_boundary', 'generate_mesh',
             'generate_case', 'run_case', 'check_result', 'gmsh', 'subprocess'}
    assert not steps & set(vars(module))


def test_paths(tmp_path, monkeypatch):                        # T8-35
    case = write(tmp_path, r0(), '算例目录')
    (tmp_path / 'elsewhere').mkdir()
    monkeypatch.chdir(tmp_path / 'elsewhere')
    for path in (case, str(case), os.path.join('..', '算例目录')):
        assert judge_result(path)['结论'] == '通过'
    assert os.listdir(tmp_path / 'elsewhere') == []


def test_command_line(tmp_path, capsys):                      # T8-36
    assert _命令行([str(write(tmp_path, r0(), 'ok'))]) == 0
    assert '通过' in capsys.readouterr().out
    assert _命令行([str(write(tmp_path, r0(['流量守恒', '压差与一维解一致']), 'redo'))]) == 1
    out = capsys.readouterr().out
    assert '流量守恒不通过，回到第 5、4 步' in out and '压差与一维解一致不通过，回到第 2、4 步' in out
    assert '加密网格' in out and '返工，回到第 2、4、5 步' in out
    (tmp_path / 'none').mkdir()
    assert _命令行([str(tmp_path / 'none')]) == 1
    captured = capsys.readouterr()
    assert captured.out == '' and '不能判定' in captured.err and '第 7 步' in captured.err
    assert _命令行([]) == 2 and _命令行(['a', 'b']) == 2 and _命令行(['--对照一维']) == 2
    env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
    for args, code in (([str(tmp_path / 'ok')], 0), ([str(tmp_path / 'redo')], 1), ([], 2)):
        done = subprocess.run([sys.executable, module.__file__] + args, env=env,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert done.returncode == code


# ---------- 可重复性 ----------

def test_repeatable(tmp_path):                                # T8-37、T8-38、T8-39
    assert judge_result(r0(项目)) == judge_result(r0(项目))
    case = write(tmp_path, r0(['流量守恒']))
    assert judge_result(case)['结论'] == '返工'
    write(tmp_path, r0())
    assert judge_result(case) == {'结论': '通过', '回到步骤': [], '对照一维': True, '原因': []}
    with pytest.raises(判定错误):
        judge_result(tmp_path / '没有')
    assert judge_result(r0())['结论'] == '通过'


# ---------- OF: 用容器真的运行 OpenFOAM ----------

def docker_ok():
    if shutil.which('docker') is None:
        return False
    return subprocess.run(['docker', 'info'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def test_real_case_with_too_few_steps_is_sent_back(tmp_path):         # T8-15
    if not docker_ok():
        if os.environ.get('CI'):
            pytest.fail('CI 里不允许跳过 OF 类测试: 不能运行 docker')
        pytest.skip('不能运行 docker，OF 类测试跳过（跳过不算通过）')
    from check_result import check_result
    from generate_boundary import generate_boundary
    from generate_case import generate_case
    from generate_mesh import generate_mesh
    from generate_structure import generate_step, 默认参数
    from run_case import run_case
    from validate_boundary import validate_boundary
    step = generate_step('S1', 默认参数('S1'), tmp_path, 'duct')
    condition = {
        '入口': {'速度': 5.0, '温度': 600.0},
        '出口': {'表压': -200.0},
        '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                 '湍流粘度': 0.002, '湍流普朗特数': 0.85},
        '壁面': '滑移',
        '多孔': {'粘性阻力系数': 1000.0, '惯性阻力系数': 20.0, '横向倍数': 10.0,
                 '换热': {'方式': '定热源', '体积热源': -50000.0}},
    }
    boundary = generate_boundary(step, condition, tmp_path)
    validate_boundary(boundary)
    mesh = generate_mesh(step, boundary, tmp_path)
    case = Path(generate_case(step, boundary, mesh, tmp_path / '算例', 步数=20)['算例目录'])
    result = check_result(run_case(case))
    old = files(case)
    verdict = judge_result(result)
    assert verdict == judge_result(case)
    assert verdict['结论'] == '返工' and 5 in verdict['回到步骤'] and verdict['对照一维'] is False
    assert '监测值稳定不通过' in [r['现象'] for r in verdict['原因']]
    assert files(case) == old
