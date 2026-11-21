# -*- coding: utf-8 -*-
'''计算任务队列: 网页提交任务，后台进程（worker.py）按提交顺序逐个计算。多人同时使用时互不干扰。

任务放在 <工作根>/任务/<任务号>/:
    任务.json        状态和元数据（号、用户、说明、状态、提交/开始/结束时间、设置、结论、说明文字）
    输入/            model.step、工况.json、标注.json（提交时复制进来，之后不再依赖网页会话）
    输出/            auto_run 的输出目录（自动运行.json、自动运行报告.txt、第NN次/...）
    日志.txt         auto_run 的逐行进展
    取消             存在时表示用户要求取消

状态: 排队 → 运行中 → 完成 / 出错 / 已取消。完成 指 auto_run 正常返回（结论可能是 通过、网格无关性未确认、需要人工、未通过）；
出错 指程序异常或计算进程意外退出。状态变化都在 <工作根>/任务/.lock 文件锁里做，并用临时文件 + 改名写入，
网页和后台进程同时读写也不会读到半个文件。

工作根: 环境变量 ZXDCFD_WORKDIR，没有时为仓库根目录下的 运行目录/。
'''
import contextlib
import fcntl
import json
import os
import shutil
import time
import uuid

仓库 = os.path.dirname(os.path.abspath(__file__))
状态集 = ('排队', '运行中', '完成', '出错', '已取消')
未结束 = ('排队', '运行中')
任务文件名 = '任务.json'
日志文件名 = '日志.txt'
取消文件名 = '取消'


class 任务错误(ValueError):
    '''提交或操作任务时的错误（找不到任务、超过限额等），信息直接给用户看。'''


def 工作根():
    return os.path.abspath(os.environ.get('ZXDCFD_WORKDIR') or os.path.join(仓库, '运行目录'))


def 任务根():
    路径 = os.path.join(工作根(), '任务')
    os.makedirs(路径, exist_ok=True)
    return 路径


def 任务目录(号):
    if not isinstance(号, str) or not 号 or not all(c.isalnum() or c == '-' for c in 号):
        raise 任务错误('任务号不合法: %r' % (号,))
    return os.path.join(任务根(), 号)


@contextlib.contextmanager
def _锁():
    with open(os.path.join(任务根(), '.lock'), 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def _写json(路径, 对象):
    临时 = '%s.%s.tmp' % (路径, uuid.uuid4().hex[:8])
    with open(临时, 'w', encoding='utf-8') as f:
        json.dump(对象, f, ensure_ascii=False, indent=1, default=str)
        f.write('\n')
    os.replace(临时, 路径)


def _读json(路径):
    with open(路径, encoding='utf-8') as f:
        return json.load(f)


def 读(号):
    路径 = os.path.join(任务目录(号), 任务文件名)
    try:
        return _读json(路径)
    except FileNotFoundError:
        raise 任务错误('找不到任务 %s（可能已过期被清理）' % 号) from None


def 改状态(号, **改):
    '''在锁里读-改-写 任务.json，返回新内容。'''
    with _锁():
        任务 = 读(号)
        任务.update(改)
        _写json(os.path.join(任务目录(号), 任务文件名), 任务)
        return 任务


def 列出(用户=None):
    '''返回任务列表，新提交的在前。用户 为 None 时列出全部。'''
    结果 = []
    for 名 in os.listdir(任务根()):
        路径 = os.path.join(任务根(), 名, 任务文件名)
        if not os.path.isfile(路径):
            continue
        try:
            任务 = _读json(路径)
        except (OSError, ValueError):
            continue
        if 用户 is None or 任务.get('用户') == 用户:
            结果.append(任务)
    return sorted(结果, key=lambda t: t['提交时间'], reverse=True)


def 排队位置(号):
    '''排队中的任务前面还有几个任务在排队（0 表示下一个就轮到它）；不在排队时返回 None。'''
    全部 = sorted(列出(), key=lambda t: t['提交时间'])
    队列 = [t['号'] for t in 全部 if t['状态'] == '排队']
    return 队列.index(号) if 号 in 队列 else None


def 提交(用户, step路径, 工况, 标注, 设置, 说明='', 每人上限=None):
    '''复制输入，建任务目录，状态为 排队。返回任务号。

    设置: auto_run 的关键字参数（网格尺寸、步数、对照一维、网格检查、运行方式、时限、单元上限、最多返工）。
    每人上限: 同一用户同时 排队 + 运行中 的任务数上限，None 表示不限。'''
    if 每人上限 is not None:
        进行中 = [t for t in 列出(用户) if t['状态'] in 未结束]
        if len(进行中) >= 每人上限:
            raise 任务错误('你已有 %d 个任务在排队或计算中（每人最多 %d 个），请等它们结束或取消后再提交'
                           % (len(进行中), 每人上限))
    号 = time.strftime('%Y%m%d-%H%M%S') + '-' + uuid.uuid4().hex[:6]
    目录 = 任务目录(号)
    os.makedirs(os.path.join(目录, '输入'))
    shutil.copyfile(step路径, os.path.join(目录, '输入', 'model.step'))
    _写json(os.path.join(目录, '输入', '工况.json'), 工况)
    _写json(os.path.join(目录, '输入', '标注.json'), 标注)
    with _锁():
        _写json(os.path.join(目录, 任务文件名),
                {'号': 号, '用户': 用户, '说明': 说明, '状态': '排队', '提交时间': time.time(),
                 '开始时间': None, '结束时间': None, '设置': 设置, '结论': None, '结论说明': '', '进程': None})
    return 号


def 领取():
    '''后台进程调用: 把最早的排队任务改成 运行中，返回任务号；没有时返回 None。'''
    with _锁():
        队列 = sorted((t for t in 列出() if t['状态'] == '排队'), key=lambda t: t['提交时间'])
        if not 队列:
            return None
        任务 = 队列[0]
        任务.update(状态='运行中', 开始时间=time.time())
        _写json(os.path.join(任务目录(任务['号']), 任务文件名), 任务)
        return 任务['号']


def 取消(号):
    '''排队中的直接取消；运行中的写取消标记，由后台进程结束计算进程后改为 已取消。'''
    with _锁():
        任务 = 读(号)
        if 任务['状态'] == '排队':
            任务.update(状态='已取消', 结束时间=time.time(), 结论说明='排队时取消')
            _写json(os.path.join(任务目录(号), 任务文件名), 任务)
        elif 任务['状态'] == '运行中':
            open(os.path.join(任务目录(号), 取消文件名), 'w').close()
        return 任务


def 要求取消(号):
    return os.path.exists(os.path.join(任务目录(号), 取消文件名))


def 日志(号, 末尾行数=None):
    try:
        with open(os.path.join(任务目录(号), 日志文件名), encoding='utf-8') as f:
            行 = f.read().splitlines()
    except FileNotFoundError:
        return ''
    return '\n'.join(行[-末尾行数:] if 末尾行数 else 行)


def 总结(号):
    '''auto_run 写的 自动运行.json；还没有时返回 None。'''
    try:
        return _读json(os.path.join(任务目录(号), '输出', '自动运行.json'))
    except (FileNotFoundError, ValueError):
        return None


def 删除(号):
    with _锁():
        if 读(号)['状态'] in 未结束:
            raise 任务错误('任务还在排队或计算中，请先取消')
        shutil.rmtree(任务目录(号), ignore_errors=True)


def 清理(保留天数):
    '''删除结束超过 保留天数 天的任务，以及 <工作根>/会话/ 下超过 1 天的临时目录。返回删掉的任务号。'''
    现在 = time.time()
    删掉 = []
    for 任务 in 列出():
        if 任务['状态'] not in 未结束 and 现在 - (任务.get('结束时间') or 任务['提交时间']) > 保留天数 * 86400:
            shutil.rmtree(任务目录(任务['号']), ignore_errors=True)
            删掉.append(任务['号'])
    会话根 = os.path.join(工作根(), '会话')
    if os.path.isdir(会话根):
        for 名 in os.listdir(会话根):
            路径 = os.path.join(会话根, 名)
            try:
                if 现在 - os.path.getmtime(路径) > 86400:
                    shutil.rmtree(路径, ignore_errors=True)
            except OSError:
                pass
    return 删掉
