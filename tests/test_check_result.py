'''docs/07-检查结果.md 的“造数”类测试点。OF 类（T7-09 至 T7-13、T7-23、T7-33、T7-35）不在这里。'''
import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from check_result import check_result, 结果错误

脚本 = Path(__file__).resolve().parent.parent / 'check_result.py'
监测量 = ['inletFlow', 'inletP', 'outletFlow', 'outletT', 'porousInP', 'porousOutP']
基本监测量 = 监测量[:4]
信息 = {'密度': 0.5, '比热': 1100.0, '入口温度': 600.0, '出口表压': -200.0, '理论入口体积流量': 10.0,
        '多孔压降估算': 150.0, '出口温度估算': 580.0, '步数': 100, '监测量': 监测量}
基准 = {'inletP': -100, 'inletFlow': -10, 'outletFlow': 10, 'outletT': 580, 'porousInP': 200, 'porousOutP': -100}
四项 = ['流量守恒', '入口流量', '监测值稳定', '残差收敛']
六项 = 四项 + ['压差与一维解一致', '出口温度与一维解一致']
方程 = ['Ux', 'Uy', 'Uz', 'p', 'T']
键 = ['通过', '迭代步数', '进出口压差', '多孔区压差', '入口流量', '出口流量', '出口温度', '吸热量',
      '多孔压降估算', '出口温度估算', '检查']
数值 = {'迭代步数': 100, '进出口压差': 150, '多孔区压差': 150, '入口流量': 10, '出口流量': 10, '出口温度': 580,
        '吸热量': 110000, '多孔压降估算': 150, '出口温度估算': 580}
无估算 = {'多孔压降估算': None, '出口温度估算': None}


def write_monitor(case, name, values):
    folder = case / 'postProcessing' / name / '0'
    folder.mkdir(parents=True, exist_ok=True)
    lines = ['# Region type : patch', '# Time        value'] + ['%d %s' % (i, v) for i, v in enumerate(values, 1)]
    (folder / 'surfaceFieldValue.dat').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def log_text(steps=100, residual=None):
    '''仿 simpleFoam 的日志: 每步各方程一行，p 解两次（第二次 0.5，不该被取到）。
    residual 是 {方程: {步下标: 值}}，值可以是文字（如 'nan'）。'''
    table = {name: {index % steps: value for index, value in (residual or {}).get(name, {}).items()}
             for name in 方程}
    lines = ['Starting time loop', '']
    for step in range(steps):
        lines.append('Time = %d' % (step + 1))
        lines.append('')
        for name in 方程:
            value = table[name].get(step, 1e-05)
            solver = 'GAMG' if name == 'p' else 'smoothSolver'
            lines.append('%s:  Solving for %s, Initial residual = %s, Final residual = 1e-07, No Iterations 2'
                         % (solver, name, value))
            if name == 'p':
                lines.append('GAMG:  Solving for p, Initial residual = 0.5, Final residual = 1e-07, No Iterations 2')
        lines.append('ExecutionTime = 1 s  ClockTime = 1 s')
        lines.append('')
    lines.append('End')
    return '\n'.join(lines) + '\n'


def make_case(root, name='case', info=None, drop=(), rows=100, values=None, edit=None, lengths=None, skip=(),
              log=True):
    '''写出基准数据 D0。info 覆盖算例信息，drop 去掉键；values 改全程的值；
    edit 是 {监测量: {行下标: 值}}，-51 就是倒数第 51 行；lengths 是 {监测量: 行数}；skip 不写数据文件。
    log 为 True 写 rows 步残差都是 1e-05 的 log.simpleFoam，为文字时照写，为 None 不写。'''
    case = Path(root) / name
    case.mkdir(parents=True)
    data = dict(copy.deepcopy(信息), **(info or {}))
    for key in drop:
        del data[key]
    (case / '算例信息.json').write_text(json.dumps(data, ensure_ascii=False), encoding='utf-8')
    for monitor in 监测量:
        if monitor in skip:
            continue
        column = [dict(基准, **(values or {}))[monitor]] * (lengths or {}).get(monitor, rows)
        for index, value in (edit or {}).get(monitor, {}).items():
            column[index] = value
        write_monitor(case, monitor, column)
    if log is not None:
        (case / 'log.simpleFoam').write_text(log_text(rows) if log is True else log, encoding='utf-8')
    return case


def failed(result):
    return [item['项目'] for item in result['检查'] if not item['通过']]


def note(result, name):
    return next(item['说明'] for item in result['检查'] if item['项目'] == name)


def saved(case):
    return json.loads((case / '结果.json').read_text(encoding='utf-8'))


def snapshot(case):
    return {str(p.relative_to(case)): p.read_bytes() for p in case.rglob('*')
            if p.is_file() and p.name != '结果.json'}


def check_numbers(result):
    for key, value in 数值.items():
        assert result[key] == pytest.approx(value, rel=1e-9), key


# ---------- 3.1 全部通过 ----------

@pytest.mark.parametrize('compare, names', [(False, 四项), (True, 六项)])
def test_d0_passes(tmp_path, compare, names):                                   # T7-01 T7-02 T7-03
    case = make_case(tmp_path)
    result = check_result(case, 对照一维=compare)
    assert [item['项目'] for item in result['检查']] == names
    assert failed(result) == [] and result['通过'] is True
    assert (case / '结果.json').is_file()
    check_numbers(result)


@pytest.mark.parametrize('change', [
    {'values': {'outletFlow': 10.009}},
    {'values': {'inletFlow': -10.09, 'outletFlow': 10.09}},
    {'edit': {'inletP': {-51: -100.09}}},
    {'edit': {'outletT': {-51: 580.009}}},
    {'edit': {'inletP': {-26: -100.09}}},
    {'edit': {'outletT': {-26: 580.009}}},
    {'edit': {'inletP': {-52: -150}, 'outletT': {-52: 590}}},
])
def test_inside_basic_thresholds(tmp_path, change):                             # T7-04
    result = check_result(make_case(tmp_path, **change))
    assert failed(result) == [] and result['通过'] is True


@pytest.mark.parametrize('values, drop', [
    ({'porousInP': 208}, 154), ({'porousInP': 192}, 146), ({'outletT': 580.5}, 150), ({'outletT': 579.5}, 150),
])
def test_inside_1d_thresholds(tmp_path, values, drop):                          # T7-05
    result = check_result(make_case(tmp_path, values=values), 对照一维=True)
    assert len(result['检查']) == 6 and failed(result) == [] and result['通过'] is True
    assert result['多孔区压差'] == pytest.approx(drop, rel=1e-9)


def test_61_rows_is_enough(tmp_path):                                           # T7-06
    result = check_result(make_case(tmp_path, rows=61))
    assert failed(result) == [] and result['迭代步数'] == 61


def test_early_transient_is_ignored(tmp_path):                                  # T7-07
    edit = {'inletP': {i: -100 * i / 49 for i in range(49)},
            'outletT': {i: 600 - 20 * i / 49 for i in range(49)}}
    result = check_result(make_case(tmp_path, edit=edit), 对照一维=True)
    assert failed(result) == [] and result['通过'] is True
    check_numbers(result)


def test_no_estimate_without_comparison(tmp_path):                              # T7-08
    case = make_case(tmp_path, info=无估算)
    result = check_result(case)
    assert failed(result) == [] and len(result['检查']) == 4
    assert saved(case)['多孔压降估算'] is None and saved(case)['出口温度估算'] is None


# ---------- 3.2 检查不通过 ----------

def check_fails(case, names, compare=False):
    result = check_result(case, 对照一维=compare)
    assert failed(result) == names and result['通过'] is False
    assert len(result['检查']) == (6 if compare else 4)
    assert saved(case)['通过'] is False
    return result


@pytest.mark.parametrize('flow', [10.011, 9.989])
def test_flow_not_conserved(tmp_path, flow):                                    # T7-14
    result = check_fails(make_case(tmp_path, values={'outletFlow': flow}), ['流量守恒'])
    assert '入口流量 10 ' in note(result, '流量守恒') and '出口流量 %s ' % flow in note(result, '流量守恒')


@pytest.mark.parametrize('flow', [10.11, 9.89])
def test_inlet_flow_off(tmp_path, flow):                                        # T7-15
    result = check_fails(make_case(tmp_path, values={'inletFlow': -flow, 'outletFlow': flow}), ['入口流量'])
    assert '计算值 %s ' % flow in note(result, '入口流量') and '理论值 10 ' in note(result, '入口流量')


def test_reversed_flow(tmp_path):                                               # T7-16
    check_fails(make_case(tmp_path, values={'inletFlow': 10, 'outletFlow': -10}), ['入口流量'])


@pytest.mark.parametrize('rows', [60, 20])
def test_too_few_rows(tmp_path, rows):                                          # T7-17
    result = check_fails(make_case(tmp_path, rows=rows), ['监测值稳定'])
    assert '共 %d 步' % rows in note(result, '监测值稳定') and result['迭代步数'] == rows


@pytest.mark.parametrize('edit, pressure, temperature', [
    ({'inletP': {-51: -100.11}}, '0.11', '0'), ({'outletT': {-51: 580.011}}, '0', '0.011'),
    ({'inletP': {-26: -100.11}}, '0.11', '0'), ({'outletT': {-26: 580.011}}, '0', '0.011'),
    ({'inletP': {-3: -99.94, -2: -100.06}}, '0.12', '0'),
    ({'outletT': {-3: 579.994, -2: 580.006}}, '0', '0.012'),
])
def test_not_steady(tmp_path, edit, pressure, temperature):                     # T7-18
    text = note(check_fails(make_case(tmp_path, edit=edit), ['监测值稳定']), '监测值稳定')
    assert '共 100 步' in text and '最后 51 行' in text
    assert 'inletP 变化 %s（' % pressure in text and 'outletT 变化 %s K' % temperature in text


@pytest.mark.parametrize('pressure, drop', [(210, 155), (190, 145)])
def test_pressure_drop_off(tmp_path, pressure, drop):                           # T7-19 T7-22
    case = make_case(tmp_path, values={'porousInP': pressure})
    result = check_fails(case, ['压差与一维解一致'], compare=True)
    assert '计算值 %d Pa' % drop in note(result, '压差与一维解一致')
    assert '估算值 150 Pa' in note(result, '压差与一维解一致')
    result = check_result(case)
    assert [item['项目'] for item in result['检查']] == 四项 and result['通过'] is True


@pytest.mark.parametrize('temperature', [580.7, 579.3, 590])
def test_outlet_temperature_off(tmp_path, temperature):                         # T7-20
    case = make_case(tmp_path, values={'outletT': temperature})
    text = note(check_fails(case, ['出口温度与一维解一致'], compare=True), '出口温度与一维解一致')
    assert '计算值 %s K' % temperature in text and '估算值 580 K' in text and '允许偏差 0.6 K' in text


def test_every_item_is_judged(tmp_path):                                        # T7-21
    result = check_fails(make_case(tmp_path, values={'outletFlow': 10.011}, rows=20), ['流量守恒', '监测值稳定'])
    assert all(item['说明'] for item in result['检查'])


@pytest.mark.parametrize('residual', [
    {'p': {-51: 9.9e-4}}, {'T': {-1: 9.9e-4}}, {'Ux': {-52: 0.5}, 'p': {-52: 0.5}, 'T': {-52: 0.5}},
])
def test_residual_inside_threshold(tmp_path, residual):                         # T7-41
    result = check_result(make_case(tmp_path, log=log_text(100, residual)))
    assert failed(result) == [] and result['通过'] is True


@pytest.mark.parametrize('residual, bad', [
    ({'p': {-51: 1.1e-3}}, 'p'), ({'T': {-26: 1.1e-3}}, 'T'), ({'Ux': {-1: 1.1e-3}}, 'Ux'),
])
def test_residual_not_converged(tmp_path, residual, bad):                       # T7-42
    text = note(check_fails(make_case(tmp_path, log=log_text(100, residual)), ['残差收敛']), '残差收敛')
    assert '日志共 100 步' in text and '超限: %s 0.0011' % bad in text


def test_residual_window_covers_short_run(tmp_path):                            # T7-43
    result = check_fails(make_case(tmp_path, rows=20, log=log_text(20, {'p': {0: 1}})), ['监测值稳定', '残差收敛'])
    assert '日志共 20 步' in note(result, '残差收敛') and '超限: p 1' in note(result, '残差收敛')


@pytest.mark.parametrize('outlet, others, passed', [
    (600.0, {}, True), (600.009, {}, True), (600.0, {'p': {-1: 1.1e-3}}, False), (599.98, {}, False),
])
def test_uniform_temperature_skips_T(tmp_path, outlet, others, passed):          # T7-45
    residual = dict({'T': {i: 0.7 for i in range(100)}}, **others)
    result = check_result(make_case(tmp_path, values={'outletT': outlet}, log=log_text(100, residual)))
    text = note(result, '残差收敛')
    assert failed(result) == ([] if passed else ['残差收敛'])
    assert ('T 不参与判断' in text) is (outlet != 599.98)
    assert ('超限: p 0.0011' in text) is bool(others) and ('超限: T 0.7' in text) is (outlet == 599.98)


def test_temperature_skip_looks_at_whole_window(tmp_path):                      # T7-45
    log = log_text(100, {'T': {i: 0.7 for i in range(100)}})
    result = check_result(make_case(tmp_path, values={'outletT': 600.011}, log=log))
    assert failed(result) == ['残差收敛'] and 'T 不参与' not in note(result, '残差收敛')
    result = check_result(make_case(tmp_path, 'b', values={'outletT': 600.0}, edit={'outletT': {-52: 590}}, log=log))
    assert failed(result) == [] and 'T 不参与判断' in note(result, '残差收敛')


@pytest.mark.parametrize('log, match', [
    (None, '没有 log.simpleFoam'),
    ('Starting time loop\nEnd\n', '没有残差行'),
    (log_text(100, {'p': {49: 'nan'}}), '第 50 个时间步 p 的残差不是有限数'),
])
def test_residual_missing(tmp_path, log, match):                                # T7-44
    case = make_case(tmp_path, log=log)
    with pytest.raises(结果错误, match=match + '.*运行第 6 步') as error:
        check_result(case)
    assert '残差' in str(error.value)
    assert '第 5 步' not in str(error.value) and '监测量' not in str(error.value)   # 第 8 步按这两个词判回第 5 步
    assert not (case / '结果.json').exists()


# ---------- 3.3 不能检查 ----------

def check_raises(case, match, compare=False):
    with pytest.raises(结果错误, match=match):
        check_result(case, 对照一维=compare)
    assert not (Path(case) / '结果.json').exists()


def test_missing_case_or_info(tmp_path):                                        # T7-24
    check_raises(tmp_path / 'none', '找不到算例目录')
    case = make_case(tmp_path)
    (case / '算例信息.json').unlink()
    check_raises(case, '没有 算例信息.json')
    (case / '算例信息.json').write_text('{密度', encoding='utf-8')
    check_raises(case, '算例信息.json 读不了')


@pytest.mark.parametrize('key', ['密度', '比热', '入口温度', '出口表压', '理论入口体积流量', '监测量'])
def test_missing_key(tmp_path, key):                                            # T7-25
    check_raises(make_case(tmp_path, drop=[key]), '没有 %s.*重新运行第 5 步' % key)


def test_missing_monitor_data(tmp_path):                                        # T7-26
    case = make_case(tmp_path, 'a')
    shutil.rmtree(case / 'postProcessing')
    check_raises(case, 'inletFlow.*先运行第 6 步')
    check_raises(make_case(tmp_path, 'b', skip=['outletT']), '没有监测量 outletT 的数据文件.*先运行第 6 步')
    check_raises(make_case(tmp_path, 'c', lengths={'inletP': 0}), '监测量 inletP 没有数据行.*先运行第 6 步')


@pytest.mark.parametrize('edit, match', [
    ({'outletT': {-1: 'nan'}}, '监测量 outletT 第 100 行'),
    ({'inletP': {49: 'inf'}}, '监测量 inletP 第 50 行'),
    ({'outletFlow': {9: 'abc'}}, '监测量 outletFlow 第 10 行'),
])
def test_value_not_finite(tmp_path, edit, match):                               # T7-27
    check_raises(make_case(tmp_path, edit=edit), match)


def test_row_counts_differ(tmp_path):                                           # T7-28
    check_raises(make_case(tmp_path, lengths={'outletFlow': 99}), '行数不一致.*outletFlow 99 行')


def test_comparison_without_estimate(tmp_path):                                 # T7-29
    check_raises(make_case(tmp_path, info=无估算), '没有一维估算值.*不对照', compare=True)


def test_inlet_outlet_drop_is_not_compared(tmp_path):                           # T7-46
    '''进出口静压差含动压变化和转弯损失，一维对照只看多孔区前后截面的压差。'''
    result = check_result(make_case(tmp_path, values={'inletP': -50}), 对照一维=True)
    assert failed(result) == [] and result['进出口压差'] == pytest.approx(175)
    assert result['多孔区压差'] == pytest.approx(150)
    assert '计算值 150 Pa' in note(result, '压差与一维解一致')


def test_comparison_without_section_monitors(tmp_path):                         # T7-47
    '''第 5 步旧版写的算例没有多孔区截面监测: 不对照时照常检查，对照时报错并提示重写算例。'''
    info = {'监测量': 基本监测量}
    result = check_result(make_case(tmp_path, 'a', info=info, skip=['porousInP', 'porousOutP']))
    assert failed(result) == [] and result['多孔区压差'] is None
    check_raises(make_case(tmp_path, 'b', info=info, skip=['porousInP', 'porousOutP']),
                 'porousInP.*重新运行第 5 步', compare=True)


@pytest.mark.parametrize('compare', ['是', 'True', 1, 0, None])
def test_invalid_comparison_argument(tmp_path, compare):                        # T7-30
    check_raises(make_case(tmp_path), '对照一维', compare=compare)
    check_raises(tmp_path / 'none', '对照一维', compare=compare)


def test_old_result_is_removed_on_error(tmp_path):                              # T7-31
    case = make_case(tmp_path, skip=['outletT'])
    (case / '结果.json').write_text('{"通过": true}', encoding='utf-8')
    check_raises(case, 'outletT')


# ---------- 3.4 输入与产物 ----------

@pytest.mark.parametrize('change, compare', [({}, True), ({'values': {'outletFlow': 10.011}, 'rows': 20}, False)])
def test_result_file_matches_return_value(tmp_path, change, compare):           # T7-32
    case = make_case(tmp_path, **change)
    result = check_result(case, 对照一维=compare)
    assert result['算例目录'] == str(case) and result['结果路径'] == str(case / '结果.json')
    text = (case / '结果.json').read_bytes().decode('utf-8')
    data = json.loads(text, parse_constant=lambda name: pytest.fail('结果.json 里有 ' + name))
    assert data == {key: value for key, value in result.items() if key not in ('算例目录', '结果路径')}
    assert list(data) == 键 and '"通过"' in text and '\\u' not in text
    assert all(sorted(item) == sorted(['项目', '通过', '说明']) for item in data['检查'])
    assert all(isinstance(item['通过'], bool) and item['说明'] for item in data['检查'])


@pytest.mark.parametrize('change, outcome', [
    ({}, '通过'), ({'values': {'outletFlow': 10.011}}, '不通过'), ({'skip': ['outletT']}, '报错'),
])
def test_inputs_are_not_modified(tmp_path, change, outcome):                    # T7-34
    case = make_case(tmp_path, **change)
    before = snapshot(case)
    argument = {'算例目录': str(case), '完成步数': 100, '监测数据': {'inletP': 'x'}}
    original = copy.deepcopy(argument)
    if outcome == '报错':
        with pytest.raises(结果错误):
            check_result(argument)
    else:
        assert check_result(argument)['通过'] is (outcome == '通过')
    assert argument == original and snapshot(case) == before
    assert (case / '结果.json').exists() is (outcome != '报错')
    assert sorted(p.name for p in case.iterdir()) == sorted(
        ['log.simpleFoam', 'postProcessing', '算例信息.json'] + (['结果.json'] if outcome != '报错' else []))


@pytest.mark.parametrize('relative', [True, False])
def test_paths_and_working_directory(tmp_path, monkeypatch, relative):          # T7-36
    case = make_case(tmp_path, '中文 算例')
    elsewhere = tmp_path / '别处'
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    result = check_result(os.path.relpath(case) if relative else str(case), 对照一维=True)
    assert result['通过'] is True and result['结果路径'] == str(case / '结果.json')
    assert (case / '结果.json').is_file() and list(elsewhere.iterdir()) == []


def run_cli(*arguments):
    return subprocess.run([sys.executable, str(脚本)] + [str(a) for a in arguments], capture_output=True,
                          encoding='utf-8', env=dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8'))


def test_command_line(tmp_path):                                                # T7-37
    done = run_cli(make_case(tmp_path, 'a'))
    assert done.returncode == 0 and done.stdout.count('通过  ') == 4 and '未通过' not in done.stdout
    done = run_cli(make_case(tmp_path, 'b'), '--对照一维')
    assert done.returncode == 0 and done.stdout.count('通过  ') == 6 and '出口温度与一维解一致' in done.stdout
    done = run_cli(make_case(tmp_path, 'c', values={'outletFlow': 10.011}))
    assert done.returncode == 1 and '未通过  流量守恒: 入口流量 10 ' in done.stdout
    assert '通过  入口流量' in done.stdout and (tmp_path / 'c' / '结果.json').is_file()
    done = run_cli(make_case(tmp_path, 'd', skip=['outletT']))
    assert done.returncode == 1 and 'outletT' in done.stderr and done.stdout == ''
    assert not (tmp_path / 'd' / '结果.json').exists()
    assert run_cli().returncode == 2
    assert run_cli(tmp_path / 'a', '--别的').returncode == 2
    assert run_cli(tmp_path / 'a', tmp_path / 'b').returncode == 2


# ---------- 3.5 可重复性 ----------

def test_repeat_gives_identical_file(tmp_path):                                 # T7-38
    case = make_case(tmp_path)
    first = check_result(case, 对照一维=True)
    content = (case / '结果.json').read_bytes()
    assert check_result(case, 对照一维=True) == first and (case / '结果.json').read_bytes() == content


def test_rerun_leaves_nothing_from_before(tmp_path):                            # T7-39
    case = make_case(tmp_path, values={'outletFlow': 10.011})
    assert check_result(case)['通过'] is False
    write_monitor(case, 'outletFlow', [10] * 100)
    assert check_result(case)['通过'] is True
    fresh = make_case(tmp_path, 'fresh')
    check_result(fresh)
    assert (case / '结果.json').read_bytes() == (fresh / '结果.json').read_bytes()


def test_error_does_not_poison_next_run(tmp_path):                              # T7-40
    check_raises(make_case(tmp_path, 'bad', skip=['outletT']), 'outletT')
    result = check_result(make_case(tmp_path, 'good'))
    assert failed(result) == [] and len(result['检查']) == 4 and result['通过'] is True
    check_numbers(result)
