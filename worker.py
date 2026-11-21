# -*- coding: utf-8 -*-
'''后台计算进程: 从任务队列（job_queue.py）里按提交顺序领取任务，每个任务单独开一个进程运行 auto_run。

用法:
    python worker.py              常驻，循环领取任务（网页 app.py 发现没有在运行时会自动启动它）
    python worker.py --跑 任务号   只跑一个任务（由常驻进程调用，也可以手动重跑）

环境变量:
    ZXDCFD_WORKDIR     工作根，与网页相同
    ZXDCFD_WORKERS     同时计算的任务数，默认 1（OpenFOAM 每个任务单核运行，可按 CPU 核数调大）
    ZXDCFD_KEEP_DAYS   任务结束后保留几天，默认 7，之后自动删除
    ZXDCFD_AUTO_RUN    测试用: 用 “模块:函数” 代替 auto_run.auto_run

同一工作根只允许一个常驻进程（用 <工作根>/任务/.worker.lock 文件锁保证），重复启动的会直接退出。
常驻进程启动时，把上次异常退出时留下的“运行中”任务标成 出错。
取消运行中的任务时，结束该任务的计算进程组、在算例目录里运行的 OpenFOAM 进程，以及该任务的 docker 容器。
'''
import fcntl
import importlib
import json
import os
import signal
import subprocess
import sys
import time
import traceback

import job_queue as 队列

轮询间隔 = 1.0
清理间隔 = 3600


def _auto_run():
    替身 = os.environ.get('ZXDCFD_AUTO_RUN')
    if 替身:
        模块, 函数 = 替身.split(':')
        return getattr(importlib.import_module(模块), 函数)
    import auto_run
    return auto_run.auto_run


def 跑(号):
    '''在当前进程里跑一个任务（状态应已是 运行中）。正常结束返回 0。'''
    目录 = 队列.任务目录(号)
    任务 = 队列.读(号)
    with open(os.path.join(目录, 队列.日志文件名), 'a', encoding='utf-8') as 日志:
        def 记(行):
            日志.write(行 + '\n')
            日志.flush()
        try:
            with open(os.path.join(目录, '输入', '工况.json'), encoding='utf-8') as f:
                工况 = json.load(f)
            with open(os.path.join(目录, '输入', '标注.json'), encoding='utf-8') as f:
                标注 = json.load(f)
            总结 = _auto_run()(os.path.join(目录, '输入', 'model.step'), 工况, os.path.join(目录, '输出'),
                               标注=标注, 记=记, **任务['设置'])
        except Exception as 错:
            记('[程序出错] %s: %s' % (type(错).__name__, 错))
            记(traceback.format_exc())
            队列.改状态(号, 状态='出错', 结束时间=time.time(), 结论说明='%s: %s' % (type(错).__name__, 错))
            return 1
    队列.改状态(号, 状态='完成', 结束时间=time.time(), 结论=总结['结论'], 结论说明=总结['说明'])
    return 0


def _结束任务进程(号, 进程):
    try:
        os.killpg(进程.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass
    # run_case 用 start_new_session 启动 Allrun，不在上面的进程组里: 按工作目录找出来结束
    目录 = os.path.realpath(队列.任务目录(号))
    for 名 in os.listdir('/proc') if os.path.isdir('/proc') else []:
        if not 名.isdigit():
            continue
        try:
            当前 = os.path.realpath(os.readlink('/proc/%s/cwd' % 名))
        except OSError:
            continue
        if 当前 == 目录 or 当前.startswith(目录 + os.sep):
            try:
                os.kill(int(名), signal.SIGKILL)
            except OSError:
                pass
    # 容器方式: 容器名前缀里带任务号（见 run_case.容器前缀）
    try:
        ids = subprocess.run(['docker', 'ps', '-q', '--filter', 'name=zxdcfd-%s-' % 号], capture_output=True,
                             text=True, timeout=30).stdout.split()
        if ids:
            subprocess.run(['docker', 'kill'] + ids, capture_output=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        pass
    进程.wait()


def _启动(号):
    环境 = dict(os.environ, ZXDCFD_CONTAINER_PREFIX='zxdcfd-%s-' % 号)
    with open(os.path.join(队列.任务目录(号), 'worker输出.txt'), 'ab') as 输出:
        进程 = subprocess.Popen([sys.executable, os.path.abspath(__file__), '--跑', 号], cwd=队列.仓库, env=环境,
                                stdin=subprocess.DEVNULL, stdout=输出, stderr=subprocess.STDOUT,
                                start_new_session=True)
    队列.改状态(号, 进程=进程.pid)
    return 进程


def 常驻(并行数=None, 保留天数=None, 一次=False):
    '''循环领取任务。一次=True 时队列空了且没有在算的任务就返回（测试用）。拿不到常驻锁时返回 False。'''
    并行数 = int(并行数 or os.environ.get('ZXDCFD_WORKERS') or 1)
    保留天数 = float(保留天数 or os.environ.get('ZXDCFD_KEEP_DAYS') or 7)
    锁文件 = open(os.path.join(队列.任务根(), '.worker.lock'), 'a')
    try:
        fcntl.flock(锁文件, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return False
    for 任务 in 队列.列出():
        if 任务['状态'] == '运行中':
            队列.改状态(任务['号'], 状态='出错', 结束时间=time.time(), 结论说明='计算服务重启，任务中断，请重新提交')
    在算 = {}
    上次清理 = 0
    while True:
        for 号, 进程 in list(在算.items()):
            if 进程.poll() is None and 队列.要求取消(号):
                _结束任务进程(号, 进程)
            if 进程.poll() is not None:
                del 在算[号]
                任务 = 队列.读(号)
                if 队列.要求取消(号):
                    队列.改状态(号, 状态='已取消', 结束时间=time.time(), 结论说明='计算中取消')
                elif 任务['状态'] == '运行中':
                    队列.改状态(号, 状态='出错', 结束时间=time.time(),
                              结论说明='计算进程意外退出（退出码 %s），详见 worker输出.txt' % 进程.returncode)
        while len(在算) < 并行数:
            号 = 队列.领取()
            if 号 is None:
                break
            在算[号] = _启动(号)
        if time.time() - 上次清理 > 清理间隔:
            队列.清理(保留天数)
            上次清理 = time.time()
        if 一次 and not 在算 and not any(t['状态'] == '排队' for t in 队列.列出()):
            锁文件.close()
            return True
        time.sleep(轮询间隔)


def 在运行():
    '''常驻进程是否在运行（看常驻锁有没有被占用）。'''
    with open(os.path.join(队列.任务根(), '.worker.lock'), 'a') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(f, fcntl.LOCK_UN)
        return False


def 确保运行():
    '''没有常驻进程时在后台启动一个（脱离当前进程，网页重启也不影响）。'''
    if 在运行():
        return False
    with open(os.path.join(队列.任务根(), 'worker.log'), 'ab') as 输出:
        subprocess.Popen([sys.executable, os.path.abspath(__file__)], cwd=队列.仓库, stdin=subprocess.DEVNULL,
                         stdout=输出, stderr=subprocess.STDOUT, start_new_session=True)
    return True


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--跑':
        sys.exit(跑(sys.argv[2]))
    if len(sys.argv) == 1:
        sys.exit(0 if 常驻() is not False else 0)
    print(__doc__)
    sys.exit(2)
