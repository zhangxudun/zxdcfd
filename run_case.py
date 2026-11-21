# -*- coding: utf-8 -*-
'''按 docs/06-求解.md 运行 OpenFOAM 算例: 第 5 步的算例目录 -> log.*、postProcessing/。

用法:
    from run_case import run_case
    结果 = run_case(算例, 时限=1800, 运行方式='容器', 镜像='opencfd/openfoam-default:2412')
    结果['完成步数']    # 另有 算例目录、日志、步数、提前收敛、监测数据、边界

    算例: 第 5 步 generate_case 的返回值，或算例目录的路径。
    时限: 秒，大于 0。超时就结束求解并报错。
    运行方式: '容器' 用 docker 启动镜像，以 root 身份执行算例目录里的 Allrun（算例目录挂载到 /case）；
              '本机' 直接执行 Allrun，要求本机装有 OpenFOAM。

命令行:
    python run_case.py 算例目录 [时限]

通过条件: 三个日志（log.gmshToFoam、log.checkMesh、log.simpleFoam）无报错，每个监测量都有数据。
    Allrun 里 checkMesh 后面有 || true，所以不只看退出码，还要读日志。FOAM Warning 不算报错。
    残差达标提前结束算通过，结果里 提前收敛 为 True，完成步数 小于 步数。
运行前先清掉上一次的结果: log.*、postProcessing/、constant/polyMesh/ 和 0 以外的时间步目录。
任何一项不通过都抛 求解错误；已经写出的日志保留，供排查。不改动第 5 步写的文件。

这一步只判断算例跑没跑完，不判断结果对不对，那是第 7 步的事。
依赖: 容器方式要 docker，本机方式要 bash 和 OpenFOAM。
'''
import glob
import json
import math
import numbers
import os
import re
import shutil
import signal
import subprocess
import uuid


class 求解错误(RuntimeError):
    '''算例目录不完整、环境不能用、求解超时、日志有报错，或监测量没有数据。'''


默认时限 = 1800
默认镜像 = 'opencfd/openfoam-default:2412'
运行方式集 = ('容器', '本机')
命令 = ('gmshToFoam', 'checkMesh', 'simpleFoam')
必需文件 = ('Allrun', 'mesh.msh', '算例信息.json')
报错关键字 = ('FOAM FATAL ERROR', 'FOAM FATAL IO ERROR', 'command not found',
              'Floating point exception', 'Segmentation fault')
监测文件名 = 'surfaceFieldValue.dat'
_时间步目录 = re.compile(r'\d+(\.\d+)?([eE][-+]?\d+)?$')


def _日志路径(目录, 名):
    return os.path.join(目录, 'log.' + 名)


def _读文本(路径):
    with open(路径, encoding='utf-8', errors='replace') as f:
        return f.read()


def _算例目录(算例):
    if isinstance(算例, dict):
        if '算例目录' not in 算例:
            raise 求解错误('算例必须是第 5 步的返回值或算例目录的路径: 字典里没有 算例目录')
        算例 = 算例['算例目录']
    try:
        return os.path.abspath(os.fspath(算例))
    except TypeError:
        raise 求解错误('算例必须是第 5 步的返回值或算例目录的路径，实际为 %r' % (算例,)) from None


def _校验时限(值):
    if isinstance(值, bool) or not isinstance(值, numbers.Real) or not math.isfinite(值) or 值 <= 0:
        raise 求解错误('时限必须是大于 0 的秒数，实际为 %r' % (值,))
    return float(值)


def _读信息(目录):
    if not os.path.isdir(目录):
        raise 求解错误('找不到算例目录: %s' % 目录)
    for 名 in 必需文件:
        if not os.path.isfile(os.path.join(目录, 名)):
            raise 求解错误('算例目录里没有 %s: %s，请先运行第 5 步' % (名, 目录))
    try:
        信息 = json.loads(_读文本(os.path.join(目录, '算例信息.json')))
    except (OSError, ValueError) as e:
        raise 求解错误('算例信息.json 读不了: %s' % e) from e
    if (not isinstance(信息, dict) or not isinstance(信息.get('监测量'), list)
            or not isinstance(信息.get('边界'), dict)):
        raise 求解错误('算例信息.json 里没有 监测量 或 边界，请重新运行第 5 步')
    return 信息


def _清旧结果(目录):
    旧 = [os.path.join(目录, 名) for 名 in os.listdir(目录)
          if 名.startswith('log.') or 名 == 'postProcessing'
          or (_时间步目录.match(名) and float(名) != 0 and os.path.isdir(os.path.join(目录, 名)))]
    旧.append(os.path.join(目录, 'constant', 'polyMesh'))
    try:
        for 路径 in 旧:
            if os.path.isdir(路径):
                shutil.rmtree(路径)
            elif os.path.exists(路径):
                os.remove(路径)
    except OSError as e:
        raise 求解错误('清不掉上一次的结果: %s' % e) from e


def _运行(目录, 时限, 运行方式, 镜像):
    '''执行 Allrun，返回 (退出码, 屏幕输出)。超时就结束进程并抛 求解错误。'''
    容器名 = None
    if 运行方式 == '容器':
        容器名 = 'zxdcfd-' + uuid.uuid4().hex
        挂载 = ['-u', '0', '-v', 目录 + ':/case', '-w', '/case']
        命令行 = (['docker', 'run', '--rm', '--name', 容器名, '--entrypoint', '/bin/bash']
                  + 挂载 + [镜像, '/case/Allrun'])
    else:
        命令行 = ['/bin/bash', os.path.join(目录, 'Allrun')]
    try:
        进程 = subprocess.Popen(命令行, cwd=目录, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, start_new_session=True)
    except OSError as e:
        raise 求解错误('启动不了 %s: %s' % (命令行[0], e)) from e
    try:
        输出 = 进程.communicate(timeout=时限)[0]
    except subprocess.TimeoutExpired:
        安静 = {'stdin': subprocess.DEVNULL, 'stdout': subprocess.DEVNULL, 'stderr': subprocess.DEVNULL}
        if 容器名:
            subprocess.run(['docker', 'kill', 容器名], **安静)
        try:
            os.killpg(进程.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
        进程.communicate()
        if 容器名:                      # 容器被强行结束时 Allrun 来不及放开文件权限，这里补做
            try:
                subprocess.run(['docker', 'run', '--rm', '--entrypoint', '/bin/chmod'] + 挂载
                               + [镜像, '-R', 'a+rwX', '/case'], timeout=120, **安静)
            except (OSError, subprocess.TimeoutExpired):
                pass
        raise 求解错误('超时: 求解超过时限 %g 秒，已结束。日志在 %s' % (时限, 目录)) from None
    return 进程.returncode, 输出.decode('utf-8', 'replace')


def _末尾(文本, 行数=6):
    行 = [l.strip() for l in 文本.splitlines() if l.strip()]
    return ' | '.join(行[-行数:]) if 行 else '（无输出）'


def _查日志(目录, 名):
    '''日志没有问题返回 None，否则返回说明。'''
    路径 = _日志路径(目录, 名)
    if not os.path.isfile(路径):
        return '没有 %s 的日志 %s' % (名, 路径)
    行 = [l.strip() for l in _读文本(路径).splitlines() if l.strip()]
    if not 行:
        return '%s 的日志是空的: %s' % (名, 路径)
    for i, l in enumerate(行):
        for 字 in 报错关键字:
            if 字 in l:
                if 字 == 'Floating point exception' and l.startswith(
                    'trapFpe: Floating point exception trapping enabled'
                ):
                    continue
                if 字 == 'command not found':
                    return ('环境里找不到 OpenFOAM，不是算例的问题。%s 的日志 %s: %s' % (名, 路径, l))
                return '%s 的日志 %s 里有报错: %s' % (名, 路径, ' | '.join(行[i:i + 6]))
    if 名 == 'checkMesh':
        if not any('Mesh OK.' in l for l in 行):
            return ('checkMesh 不通过，日志 %s: %s。请回到第 4 步减小网格尺寸，或修几何'
                    % (路径, ' | '.join(行[-6:])))
    elif 行[-1] != 'End':
        return '%s 的日志 %s 不以 End 结束，最后几行: %s' % (名, 路径, ' | '.join(行[-6:]))
    return None


def _读边界文件(目录):
    '''读 constant/polyMesh/boundary，返回 {边界名: 类型}；文件不存在返回 None。'''
    路径 = os.path.join(目录, 'constant', 'polyMesh', 'boundary')
    if not os.path.isfile(路径):
        return None
    文本 = re.sub(r'//[^\n]*', '', re.sub(r'/\*.*?\*/', '', _读文本(路径), flags=re.S))
    边界 = {}
    for 名, 内容 in re.findall(r'(\w+)\s*\{([^{}]*)\}', 文本):
        类 = re.search(r'(?<![A-Za-z])type\s+(\w+)\s*;', 内容)
        if 名 != 'FoamFile' and 类:
            边界[名] = 类.group(1)
    return 边界


def _读监测(目录, 名, 完成步数):
    文件 = sorted(glob.glob(os.path.join(glob.escape(目录), 'postProcessing', glob.escape(名),
                                         '*', 监测文件名)))
    if not 文件:
        raise 求解错误('没有监测量 %s 的数据文件 postProcessing/%s/0/%s' % (名, 名, 监测文件名))
    行 = [l.split() for l in _读文本(文件[-1]).splitlines()
          if l.strip() and not l.lstrip().startswith('#')]
    if not 行:
        raise 求解错误('监测量 %s 没有数据行: %s' % (名, 文件[-1]))
    for i, 段 in enumerate(行, 1):
        try:
            值 = float(段[-1])
        except ValueError:
            值 = float('nan')
        if not math.isfinite(值):
            raise 求解错误('监测量 %s 第 %d 行的值不是有限数（%s）: %s' % (名, i, 段[-1], 文件[-1]))
    if len(行) != 完成步数:
        raise 求解错误('监测量 %s 有 %d 行数据，但完成步数是 %d: %s' % (名, len(行), 完成步数, 文件[-1]))
    return 文件[-1]


def run_case(算例, 时限=默认时限, 运行方式='容器', 镜像=默认镜像):
    '''运行算例并检查日志和监测数据，返回:

        {'算例目录': 绝对路径, '日志': {命令名: 路径}, '步数': 第 5 步设的步数,
         '完成步数': n, '提前收敛': 布尔, '监测数据': {监测量: 数据文件路径},
         '边界': {边界名: 类型}}      # 来自 constant/polyMesh/boundary

    不通过时抛 求解错误。不改动第 5 步写的文件。
    '''
    目录 = _算例目录(算例)
    时限 = _校验时限(时限)
    if 运行方式 not in 运行方式集:
        raise 求解错误('运行方式只能是: %s，实际为 %r' % (' / '.join(运行方式集), 运行方式))
    信息 = _读信息(目录)
    if 运行方式 == '容器' and shutil.which('docker') is None:
        raise 求解错误('找不到 docker: 容器方式要先装好 docker，或改用 运行方式=\'本机\'')
    _清旧结果(目录)

    退出码, 输出 = _运行(目录, 时限, 运行方式, 镜像)
    if 退出码 != 0:
        有 = [名 for 名 in 命令 if os.path.isfile(_日志路径(目录, 名))]
        if not 有:
            raise 求解错误('Allrun 没有跑到 gmshToFoam，退出码 %d: %s' % (退出码, _末尾(输出)))
        停 = 有[-1]
        说明 = _查日志(目录, 停) or '%s 的日志 %s 最后几行: %s' % (
            停, _日志路径(目录, 停), _末尾(_读文本(_日志路径(目录, 停))))
        raise 求解错误('停在 %s，退出码 %d。%s' % (停, 退出码, 说明))
    for 名 in 命令:
        说明 = _查日志(目录, 名)
        if 说明:
            raise 求解错误(说明)

    边界 = _读边界文件(目录)
    if 边界 is not None and 边界 != 信息['边界']:
        raise 求解错误('转出的网格边界与算例信息不一致: 网格里是 %s，算例信息里是 %s'
                       % (sorted(边界.items()), sorted(信息['边界'].items())))
    日志 = _读文本(_日志路径(目录, 'simpleFoam'))
    # 每个时间步一行 Time = …；温差换热时 deltaT 可能小于 1（见 generate_case.伪时间步），所以数行数，不读时间值
    时间 = re.findall(r'^Time = ([0-9.eE+-]+)s?\s*$', 日志, re.M)
    if not 时间:
        raise 求解错误('simpleFoam 的日志 %s 里没有时间步' % _日志路径(目录, 'simpleFoam'))
    完成步数 = len(时间)
    提前收敛 = bool(re.search(r'SIMPLE solution converged in [0-9.eE+-]+ iterations', 日志))
    步数 = 信息.get('步数')
    if 步数 is not None and not 提前收敛 and 完成步数 != 步数:
        raise 求解错误('simpleFoam 只跑了 %d 步，应为 %d 步，日志里也没有提前收敛的记录: %s'
                       % (完成步数, 步数, _日志路径(目录, 'simpleFoam')))
    监测数据 = {名: _读监测(目录, 名, 完成步数) for 名 in 信息['监测量']}
    return {'算例目录': 目录, '日志': {名: _日志路径(目录, 名) for 名 in 命令}, '步数': 步数,
            '完成步数': 完成步数, '提前收敛': 提前收敛, '监测数据': 监测数据, '边界': 边界}


if __name__ == '__main__':
    import sys
    try:
        输入时限 = float(sys.argv[2]) if len(sys.argv) == 3 else 默认时限
    except ValueError:
        输入时限 = None
    if len(sys.argv) not in (2, 3) or 输入时限 is None:
        print(__doc__)
        sys.exit(2)
    try:
        结果 = run_case(sys.argv[1], 输入时限)
    except 求解错误 as e:
        print('求解不通过: %s' % e, file=sys.stderr)
        sys.exit(1)
    print('完成 %d 步（设定 %s 步）%s，监测量 %s'
          % (结果['完成步数'], 结果['步数'], '，提前收敛' if 结果['提前收敛'] else '',
             ' / '.join(结果['监测数据'])))
    print(结果['算例目录'])
