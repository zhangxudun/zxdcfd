# -*- coding: utf-8 -*-
'''任务队列（job_queue.py）和后台计算进程（worker.py）的测试。auto_run 用 tests/stub_auto_run.py 代替，不需要 OpenFOAM。

JQ-01  提交 → 排队 → 常驻进程领取并在子进程里计算 → 完成，日志和总结都能读到，输入原样传给 auto_run。
JQ-02  按提交顺序计算；排队位置正确。
JQ-03  每人同时进行中的任务数有上限，别人的任务不算。
JQ-04  排队中取消: 直接 已取消，不会被计算。
JQ-05  计算中取消: 结束计算进程，状态 已取消。
JQ-06  auto_run 抛异常: 状态 出错，原因写进 结论说明 和日志。
JQ-07  上次异常退出留下的 运行中 任务，常驻进程启动时标成 出错。
JQ-08  同一工作根只能有一个常驻进程。
JQ-09  清理: 结束超过保留天数的任务被删除，未结束的不删。
JQ-10  任务号不合法时报错（防止路径穿越）。
'''
import json
import os
import subprocess
import sys
import time

import pytest

仓库 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, 仓库)
import job_queue as 队列  # noqa: E402
import worker  # noqa: E402


@pytest.fixture
def 环境(tmp_path, monkeypatch):
    monkeypatch.setenv('ZXDCFD_WORKDIR', str(tmp_path / '工作根'))
    monkeypatch.setenv('ZXDCFD_AUTO_RUN', 'stub_auto_run:假运行')
    monkeypatch.setenv('PYTHONPATH', os.path.join(仓库, 'tests'))
    monkeypatch.delenv('ZXDCFD_STUB_SLEEP', raising=False)
    monkeypatch.delenv('ZXDCFD_STUB_FAIL', raising=False)
    monkeypatch.setattr(worker, '轮询间隔', 0.2)
    step = tmp_path / 'model.step'
    step.write_text('不是真的 STEP，替身不读它')
    return str(step)


标注 = {'面': {'5': 'inlet', '18': 'outlet'}, '多孔体': [2]}


def 提交(step, 用户='甲', **改):
    return 队列.提交(用户, step, {'入口': {'速度': 5.0}}, 标注, dict({'步数': 900}, **改), 说明='测试')


def test_JQ_01_提交到完成(环境):
    号 = 提交(环境)
    assert 队列.读(号)['状态'] == '排队'
    assert worker.常驻(一次=True) is True
    任务 = 队列.读(号)
    assert 任务['状态'] == '完成' and 任务['结论'] == '通过' and 任务['结论说明'] == '替身说明'
    assert 任务['开始时间'] and 任务['结束时间'] >= 任务['开始时间']
    总结 = 队列.总结(号)
    assert 总结['设置'] == {'步数': 900} and 总结['工况'] == {'入口': {'速度': 5.0}}
    assert '"5": "inlet"' in 队列.日志(号) and '[结论] 通过' in 队列.日志(号)
    assert os.path.isfile(os.path.join(队列.任务目录(号), '输入', 'model.step'))


def test_JQ_02_按顺序(环境):
    号们 = [提交(环境, 用户=u) for u in ('甲', '乙', '丙')]
    assert [队列.排队位置(h) for h in 号们] == [0, 1, 2]
    worker.常驻(一次=True)
    任务们 = [队列.读(h) for h in 号们]
    assert all(t['状态'] == '完成' for t in 任务们)
    开始 = [t['开始时间'] for t in 任务们]
    assert 开始 == sorted(开始)
    assert [t['号'] for t in 队列.列出('乙')] == [号们[1]]


def test_JQ_03_每人上限(环境):
    提交(环境)
    队列.提交('甲', 环境, {}, 标注, {}, 每人上限=2)
    with pytest.raises(队列.任务错误, match='每人最多 2 个'):
        队列.提交('甲', 环境, {}, 标注, {}, 每人上限=2)
    队列.提交('乙', 环境, {}, 标注, {}, 每人上限=2)


def test_JQ_04_排队中取消(环境):
    号 = 提交(环境)
    队列.取消(号)
    assert 队列.读(号)['状态'] == '已取消'
    worker.常驻(一次=True)
    assert 队列.读(号)['状态'] == '已取消' and 队列.日志(号) == ''


def test_JQ_05_计算中取消(环境, monkeypatch):
    monkeypatch.setenv('ZXDCFD_STUB_SLEEP', '30')
    号 = 提交(环境)
    进程 = subprocess.Popen([sys.executable, '-c', 'import worker; worker.轮询间隔 = 0.2; worker.常驻(一次=True)'],
                            cwd=仓库)
    try:
        for _ in range(100):
            if 队列.读(号)['状态'] == '运行中' and 队列.日志(号):
                break
            time.sleep(0.2)
        assert 队列.读(号)['状态'] == '运行中'
        开始 = time.time()
        队列.取消(号)
        进程.wait(timeout=20)
        assert time.time() - 开始 < 20
    finally:
        进程.kill()
    任务 = 队列.读(号)
    assert 任务['状态'] == '已取消' and 任务['结论说明'] == '计算中取消'


def test_JQ_06_出错(环境, monkeypatch):
    monkeypatch.setenv('ZXDCFD_STUB_FAIL', '1')
    号 = 提交(环境)
    worker.常驻(一次=True)
    任务 = 队列.读(号)
    assert 任务['状态'] == '出错' and '替身按要求出错' in 任务['结论说明']
    assert '[程序出错]' in 队列.日志(号)


def test_JQ_07_重启时标出错(环境):
    号 = 提交(环境)
    assert 队列.领取() == 号
    worker.常驻(一次=True)
    任务 = 队列.读(号)
    assert 任务['状态'] == '出错' and '重启' in 任务['结论说明']


def test_JQ_08_只能一个常驻进程(环境, monkeypatch):
    monkeypatch.setenv('ZXDCFD_STUB_SLEEP', '2')
    提交(环境)
    进程 = subprocess.Popen([sys.executable, '-c', 'import worker; worker.常驻(一次=True)'], cwd=仓库)
    try:
        for _ in range(50):
            if worker.在运行():
                break
            time.sleep(0.1)
        assert worker.在运行()
        assert worker.常驻(一次=True) is False
        assert worker.确保运行() is False
    finally:
        进程.wait(timeout=60)
    assert not worker.在运行()


def test_JQ_09_清理(环境):
    旧, 新, 排 = 提交(环境), 提交(环境, 用户='乙'), None
    worker.常驻(一次=True)
    排 = 提交(环境, 用户='丙')
    队列.改状态(旧, 结束时间=time.time() - 8 * 86400)
    队列.改状态(排, 提交时间=time.time() - 30 * 86400)
    assert 队列.清理(7) == [旧]
    assert {t['号'] for t in 队列.列出()} == {新, 排}


def test_JQ_10_任务号不合法(环境):
    for 坏 in ('../x', '', 'a/b', None):
        with pytest.raises(队列.任务错误):
            队列.任务目录(坏)
    with pytest.raises(队列.任务错误, match='找不到'):
        队列.读('20990101-000000-abcdef')
