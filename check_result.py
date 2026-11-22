# -*- coding: utf-8 -*-
'''按 docs/07-检查结果.md 检查 OpenFOAM 算例的结果: 算例信息.json、postProcessing/、log.simpleFoam -> 结果.json。

用法:
    from check_result import check_result
    结果 = check_result(算例, 对照一维=False)
    结果['通过']    # 另有 迭代步数、进出口压差、多孔区压差、入口流量、出口流量、出口温度、吸热量、
                    # 多孔压降估算、出口温度估算、检查、算例目录、结果路径

    算例: 第 6 步 run_case 的返回值，或算例目录的路径。
    对照一维: True 时在四项基本检查之外，再拿压差和出口温度与第 5 步的一维估算对照。
              几何没有一维估算值时只能取 False。

命令行:
    python check_result.py 算例目录 [--对照一维]
    全部通过返回 0，有检查不通过或不能检查返回 1，用法不对返回 2。

检查项（阈值都是严格小于；末值是监测数据最后一行，变化量是最后 51 行里的最大值减最小值）:
    流量守恒            |inletFlow 末值 + outletFlow 末值| < 0.1% × |inletFlow 末值|
    入口流量            |入口流量 − 理论入口体积流量| < 1% × 理论入口体积流量
    监测值稳定          数据多于 60 行，inletP 变化量 < 0.1% × |inletP 末值|，outletT 变化量 < 0.01 K
    残差收敛            log.simpleFoam 里每个方程最后 51 步的初始残差最大值 < 1e-3（不足 51 步取全部）；
                        outletT 最后 51 行都等于入口温度（没有换热）时 T 不参与
    压差与一维解一致    |多孔区压差 − 多孔压降估算| < 3% × 多孔压降估算
                        多孔区压差 = (porousInP 末值 − porousOutP 末值) × 密度，取多孔区前后两个截面，
                        不含进出口动压变化和转弯、变截面的损失；算例没有这两个监测量时不能对照一维
    出口温度与一维解一致 |outletT 末值 − 出口温度估算| < 3% × |入口温度 − 出口温度估算|（相对温降）

监测值稳定看的是最后 51 行整个窗口，不是只比首尾两行: 首尾相等但中间振荡的数据判为不稳定。
残差取每个时间步里各方程第一次求解的初始残差（p 一步解几次的只取第一次）。

结局分三种:
    全部通过: 写 结果.json，通过 为 True。
    检查不通过: 数据齐全但有检查项没过。不报错，照样写 结果.json，通过 为 False，每一项都有结论。
    不能检查: 缺文件、缺数据或数据不是有限数（含没有 log.simpleFoam、日志里没有残差、残差不是有限数）。
              抛 结果错误，不写 结果.json。
运行一开始就删掉旧的 结果.json，所以报错以后目录里没有 结果.json。
只读 算例信息.json、postProcessing/ 和 log.simpleFoam，只写 结果.json。不运行 OpenFOAM。
读日志只为取残差，日志有没有报错是第 6 步的事。

这一步只判断结果对不对，不判断不通过以后回到哪一步，那是第 8 步的事。
'''
import glob
import json
import math
import numbers
import os
import re


class 结果错误(RuntimeError):
    '''算例目录或算例信息不完整、监测数据缺失或不是有限数，或要求对照一维但没有估算值。'''


结果文件名 = '结果.json'
信息文件名 = '算例信息.json'
监测文件名 = 'surfaceFieldValue.dat'
必需键 = ('密度', '比热', '入口温度', '出口表压', '理论入口体积流量')
必需监测量 = ('inletFlow', 'inletP', 'outletFlow', 'outletT')
估算键 = ('多孔压降估算', '出口温度估算')
截面监测量 = ('porousInP', 'porousOutP')
稳定间隔 = 50          # 判稳的窗口是最后 稳定间隔 + 1 行
稳定最少行数 = 60      # 数据要多于这么多行
流量守恒容差 = 1e-3
入口流量容差 = 1e-2
压力稳定容差 = 1e-3
温度稳定容差 = 0.01    # K
日志文件名 = 'log.simpleFoam'
残差窗口 = 稳定间隔 + 1   # 与“监测值稳定”同一个窗口
残差限 = 1e-3
温度方程 = 'T'
_残差行 = re.compile(r'Solving for (\w+), Initial residual = ([^,\s]+)')
一维容差 = 0.03


def _读文本(路径):
    with open(路径, encoding='utf-8', errors='replace') as f:
        return f.read()


def _算例目录(算例):
    if isinstance(算例, dict):
        if '算例目录' not in 算例:
            raise 结果错误('算例必须是第 6 步的返回值或算例目录的路径: 字典里没有 算例目录')
        算例 = 算例['算例目录']
    try:
        return os.path.abspath(os.fspath(算例))
    except TypeError:
        raise 结果错误('算例必须是第 6 步的返回值或算例目录的路径，实际为 %r' % (算例,)) from None


def _清旧结果(目录):
    路径 = os.path.join(目录, 结果文件名)
    try:
        if os.path.isfile(路径) or os.path.islink(路径):
            os.remove(路径)
    except OSError as e:
        raise 结果错误('清不掉上一次的 %s: %s' % (结果文件名, e)) from e


def _是有限数(值):
    return (not isinstance(值, bool) and isinstance(值, numbers.Real) and math.isfinite(值))


def _读信息(目录):
    if not os.path.isdir(目录):
        raise 结果错误('找不到算例目录: %s' % 目录)
    路径 = os.path.join(目录, 信息文件名)
    if not os.path.isfile(路径):
        raise 结果错误('算例目录里没有 %s: %s，请先运行第 5 步' % (信息文件名, 目录))
    try:
        信息 = json.loads(_读文本(路径))
    except (OSError, ValueError) as e:
        raise 结果错误('%s 读不了: %s' % (信息文件名, e)) from e
    if not isinstance(信息, dict):
        raise 结果错误('%s 读不了: 内容不是字典，请重新运行第 5 步' % 信息文件名)
    for 键 in 必需键 + ('监测量',):
        if 键 not in 信息:
            raise 结果错误('%s 里没有 %s，请重新运行第 5 步' % (信息文件名, 键))
    for 键 in 必需键:
        if not _是有限数(信息[键]):
            raise 结果错误('%s 里的 %s 不是有限数（%r），请重新运行第 5 步' % (信息文件名, 键, 信息[键]))
    监测量 = 信息['监测量']
    if not isinstance(监测量, list) or not all(isinstance(名, str) and 名 for 名 in 监测量):
        raise 结果错误('%s 里的 监测量 不是名字列表，请重新运行第 5 步' % 信息文件名)
    for 名 in 必需监测量:
        if 名 not in 监测量:
            raise 结果错误('%s 的 监测量 里没有 %s，请重新运行第 5 步' % (信息文件名, 名))
    for 键 in 估算键:
        if 信息.get(键) is not None and not _是有限数(信息[键]):
            raise 结果错误('%s 里的 %s 不是有限数（%r），请重新运行第 5 步' % (信息文件名, 键, 信息[键]))
    return 信息


def _读监测(目录, 名):
    '''返回监测量每一步的值。读法与第 6 步相同: # 开头的是注释，其余每行一步，值取最后一列。'''
    文件 = sorted(glob.glob(os.path.join(glob.escape(目录), 'postProcessing', glob.escape(名),
                                         '*', 监测文件名)))
    if not 文件:
        raise 结果错误('没有监测量 %s 的数据文件 postProcessing/%s/0/%s，请先运行第 6 步'
                       % (名, 名, 监测文件名))
    try:
        文本 = _读文本(文件[-1])
    except OSError as e:
        raise 结果错误('监测量 %s 的数据文件读不了: %s' % (名, e)) from e
    行 = [l.split() for l in 文本.splitlines() if l.strip() and not l.lstrip().startswith('#')]
    if not 行:
        raise 结果错误('监测量 %s 没有数据行: %s，请先运行第 6 步' % (名, 文件[-1]))
    值列 = []
    for i, 段 in enumerate(行, 1):
        try:
            值 = float(段[-1])
        except ValueError:
            值 = float('nan')
        if not math.isfinite(值):
            raise 结果错误('监测量 %s 第 %d 行数据的值不是有限数（%s）: %s' % (名, i, 段[-1], 文件[-1]))
        值列.append(值)
    return 值列


def _读残差(目录):
    '''返回 {方程: [每个时间步第一次求解的初始残差]}，方程按日志里出现的先后排。'''
    路径 = os.path.join(目录, 日志文件名)
    if not os.path.isfile(路径):
        raise 结果错误('算例目录里没有 %s，读不到残差，请先运行第 6 步' % 日志文件名)
    try:
        文本 = _读文本(路径)
    except OSError as e:
        raise 结果错误('%s 读不了，读不到残差: %s，请重新运行第 6 步' % (日志文件名, e)) from e
    残差, 本步, 步号 = {}, set(), 0
    for 行 in 文本.splitlines():
        if 行.startswith('Time = '):
            本步, 步号 = set(), 步号 + 1
            continue
        找到 = _残差行.search(行)
        if not 找到 or 找到.group(1) in 本步:
            continue
        名, 文字 = 找到.groups()
        try:
            值 = float(文字)
        except ValueError:
            值 = float('nan')
        if not math.isfinite(值):
            raise 结果错误('%s 第 %d 个时间步 %s 的残差不是有限数（%s），请重新运行第 6 步'
                           % (日志文件名, 步号, 名, 文字))
        本步.add(名)
        残差.setdefault(名, []).append(值)
    if not 残差:
        raise 结果错误('%s 里没有残差行（Solving for ...），请重新运行第 6 步' % 日志文件名)
    return 残差


def _检查残差(残差, 出口温度窗口, 入口温度):
    步数 = max(len(值列) for 值列 in 残差.values())
    窗口最大 = {名: max(值列[-残差窗口:]) for 名, 值列 in 残差.items()}
    说明 = '日志共 %d 步' % 步数
    # 没有换热时温度场处处等于入口温度，OpenFOAM 的归一化因子接近 0，T 的初始残差只是舍入噪声（常在 0.1–1），
    # 不说明没收敛。出口温度在整个窗口里都等于入口温度，就认为没有换热，T 不参与判断。
    if 温度方程 in 窗口最大 and max(abs(t - 入口温度) for t in 出口温度窗口) <= 温度稳定容差:
        del 窗口最大[温度方程]
        说明 += '；出口温度在最后 %d 行里一直等于入口温度 %s K，没有换热，%s 不参与判断' % (
            len(出口温度窗口), _数(入口温度), 温度方程)
    超限 = [名 for 名, 值 in 窗口最大.items() if not 值 < 残差限]
    说明 += ('；最后 %d 步里各方程初始残差的最大值: %s（都要小于 %s）'
           % (min(步数, 残差窗口), '，'.join('%s %s' % (名, _数(值)) for 名, 值 in 窗口最大.items()) or '无',
              _数(残差限)))
    if 超限:
        说明 += '；超限: %s' % '，'.join('%s %s' % (名, _数(窗口最大[名])) for 名 in 超限)
    return _项('残差收敛', not 超限, 说明)


热平衡容差 = 1e-2       # 烟气放热与水吸热之差 / 水吸热


def _读水侧(目录, 水侧):
    '''读 postProcessing/waterSide.dat（每步: 时间 水出口温度 水吸热 烟气按旧水温算的放热）和 waterProfile.dat。'''
    行列 = []
    for 键 in ('监测文件', '水温分布文件'):
        路径 = os.path.join(目录, *水侧[键].split('/'))
        if not os.path.isfile(路径):
            raise 结果错误('水侧耦合算例没有 %s，请先运行第 6 步' % 水侧[键])
        try:
            行 = [l.split() for l in _读文本(路径).splitlines() if l.strip() and not l.lstrip().startswith('#')]
            值 = [[float(x) for x in 段] for 段 in 行]
        except (OSError, ValueError) as e:
            raise 结果错误('%s 读不了: %s' % (水侧[键], e)) from e
        if not 值 or not all(len(段) == len(值[0]) and all(math.isfinite(x) for x in 段) for 段 in 值):
            raise 结果错误('%s 没有数据，或数据不是有限数，请重新运行第 6 步' % 水侧[键])
        行列.append(值)
    return 行列


def _检查水侧(目录, 信息, 吸热量, 对照一维):
    '''水侧耦合时的三项检查（加对照一维时的一项），返回 (检查项列表, 结果字段)。'''
    水侧 = 信息['水侧']
    记录, 分布 = _读水侧(目录, 水侧)
    窗口 = 稳定间隔 + 1
    水出口 = [r[1] for r in 记录]
    水吸热, 烟放热 = 记录[-1][2], 记录[-1][3]
    温度变化 = _窗口变化(水出口, 窗口)
    差 = abs(水吸热 - 烟放热)
    差限 = 1e-3 * abs(水吸热) + 1e-9
    检查 = [_项('水侧稳定', len(记录) >= 窗口 and 温度变化 < 温度稳定容差 and 差 < 差限,
                 '共 %d 行；最后 %d 行里水出口温度变化 %s K（允许 %s K）；最后一步水吸热 %s W，'
                 '按上一步水温算的烟气放热 %s W，相差 %s W（允许 %s W，即 0.1%%）'
                 % (len(记录), 窗口, _数(温度变化), _数(温度稳定容差), _数(水吸热), _数(烟放热), _数(差), _数(差限)))]
    绝热 = 信息.get('壁面换热', {'方式': '绝热'}).get('方式') == '绝热'
    平衡差 = abs(吸热量 - 水吸热)
    平衡限 = 热平衡容差 * abs(水吸热)
    if 绝热:
        检查.append(_项('热平衡', 平衡差 < 平衡限,
                        '烟气放热（按进出口温度）%s W，水吸热 %s W，相差 %s W（允许 %s W，即水吸热的 1%%）'
                        % (_数(吸热量), _数(水吸热), _数(平衡差), _数(平衡限))))
    else:
        检查.append(_项('热平衡', True, '外壁不是绝热，烟气还和外壁换热，不对照（烟气放热 %s W，水吸热 %s W）'
                        % (_数(吸热量), _数(水吸热))))
    最高 = max(r[1] for r in 分布)
    饱和 = 水侧.get('水饱和温度')
    if 饱和 is None:
        检查.append(_项('水未沸腾', True, '没给水饱和温度，不检查；水温最高 %s K' % _数(最高)))
    else:
        检查.append(_项('水未沸腾', 最高 < 饱和, '水温最高 %s K，饱和温度 %s K（要低于饱和温度；单相模型）'
                        % (_数(最高), _数(饱和))))
    if 对照一维:
        估算 = 水侧.get('水出口温度估算')
        if not _是有限数(估算):
            raise 结果错误('%s 的 水侧 里没有 水出口温度估算，不能对照一维，请改为不对照' % 信息文件名)
        偏差 = abs(水出口[-1] - 估算)
        限 = 一维容差 * abs(估算 - 水侧['水入口温度'])
        检查.append(_项('水出口温度与一维解一致', 偏差 < 限,
                        '计算值 %s K，估算值 %s K，相差 %s K（允许偏差 %s K，即估算温升 %s K 的 3%%）'
                        % (_数(水出口[-1]), _数(估算), _数(偏差), _数(限), _数(abs(估算 - 水侧['水入口温度'])))))
    return 检查, {'水出口温度': 水出口[-1], '水吸热量': 水吸热, '水温最高': 最高,
                  '水出口温度估算': 水侧.get('水出口温度估算')}


def _项(项目, 通过, 说明):
    return {'项目': 项目, '通过': bool(通过), '说明': 说明}


def _数(值):
    return '%.6g' % 值


def _窗口变化(值列, 窗口):
    '''最后 窗口 行里的最大值减最小值。'''
    末段 = 值列[-窗口:]
    return max(末段) - min(末段)


def _检查稳定(行数, 压力, 温度):
    窗口 = 稳定间隔 + 1
    if 行数 < 窗口:
        return _项('监测值稳定', False,
                   '共 %d 步，不足 %d 步，凑不齐最后 %d 行，压力变化量和温度变化量无法计算；数据要多于 %d 行'
                   % (行数, 窗口, 窗口, 稳定最少行数))
    压力变化 = _窗口变化(压力, 窗口)
    压力限 = 压力稳定容差 * abs(压力[-1])
    温度变化 = _窗口变化(温度, 窗口)
    通过 = 行数 > 稳定最少行数 and 压力变化 < 压力限 and 温度变化 < 温度稳定容差
    说明 = ('共 %d 步（要多于 %d 步）；最后 %d 行里最大值减最小值: inletP 变化 %s（允许 %s，即末值 %s 的 0.1%%），'
            'outletT 变化 %s K（允许 %s K）'
            % (行数, 稳定最少行数, 窗口, _数(压力变化), _数(压力限), _数(压力[-1]),
               _数(温度变化), _数(温度稳定容差)))
    return _项('监测值稳定', 通过, 说明)


def check_result(算例, 对照一维=False):
    '''检查算例结果并写 结果.json，返回:

        {'通过': 布尔, '迭代步数': n, '进出口压差': Pa, '多孔区压差': Pa 或 None, '入口流量': m³/s, '出口流量': m³/s,
         '出口温度': K, '吸热量': W, '多孔压降估算': Pa 或 None, '出口温度估算': K 或 None,
         '检查': [{'项目': 名, '通过': 布尔, '说明': 文字}, ...],
         '算例目录': 绝对路径, '结果路径': 结果.json 的路径}      # 最后两项不写进 结果.json

    检查不通过不报错，通过 为 False。不能检查时抛 结果错误，不写 结果.json。
    '''
    目录 = _算例目录(算例)
    _清旧结果(目录)
    if not isinstance(对照一维, bool):
        raise 结果错误('对照一维 只能是 True 或 False，实际为 %r' % (对照一维,))
    信息 = _读信息(目录)
    压降估算, 温度估算 = (信息.get(键) for 键 in 估算键)
    if 对照一维 and (压降估算 is None or 温度估算 is None):
        raise 结果错误('这个几何没有一维估算值（%s 里的 %s 为空），不能对照一维，请改为不对照'
                       % (信息文件名, '、'.join(键 for 键 in 估算键 if 信息.get(键) is None)))
    有截面 = all(名 in 信息['监测量'] for 名 in 截面监测量)
    if 对照一维 and not 有截面:
        raise 结果错误('%s 的 监测量 里没有多孔区前后截面的压力（%s），不能对照一维压差，请重新运行第 5 步'
                       % (信息文件名, '、'.join(截面监测量)))

    数据 = {名: _读监测(目录, 名) for 名 in 信息['监测量']}
    行数集 = {名: len(值列) for 名, 值列 in 数据.items()}
    if len(set(行数集.values())) != 1:
        raise 结果错误('监测量的行数不一致: %s，请重新运行第 6 步'
                       % '，'.join('%s %d 行' % (名, n) for 名, n in 行数集.items()))
    行数 = 行数集['inletP']
    残差 = _读残差(目录)

    密度, 比热, 入口温度, 出口表压, 理论流量 = (float(信息[键]) for 键 in 必需键)
    入口监测流量 = 数据['inletFlow'][-1]
    出口监测流量 = 数据['outletFlow'][-1]
    入口流量 = -入口监测流量 + 0.0
    出口流量 = 出口监测流量 + 0.0
    出口温度 = 数据['outletT'][-1]
    进出口压差 = 数据['inletP'][-1] * 密度 - 出口表压
    多孔区压差 = (数据['porousInP'][-1] - 数据['porousOutP'][-1]) * 密度 if 有截面 else None
    吸热量 = 密度 * 入口流量 * 比热 * (入口温度 - 出口温度)
    for 名, 值 in (('进出口压差', 进出口压差), ('多孔区压差', 多孔区压差), ('吸热量', 吸热量)):
        if 值 is None:
            continue
        if not math.isfinite(值):
            raise 结果错误('算出的 %s 不是有限数（%r），请检查 %s 和监测数据' % (名, 值, 信息文件名))

    流量差 = abs(入口监测流量 + 出口监测流量)
    流量限 = 流量守恒容差 * abs(入口监测流量)
    入口偏差 = abs(入口流量 - 理论流量)
    入口限 = 入口流量容差 * 理论流量
    检查 = [
        _项('流量守恒', 流量差 < 流量限,
            '入口流量 %s m³/s，出口流量 %s m³/s，相差 %s（允许 %s，即入口流量的 0.1%%）'
            % (_数(入口流量), _数(出口流量), _数(流量差), _数(流量限))),
        _项('入口流量', 入口偏差 < 入口限,
            '计算值 %s m³/s，理论值 %s m³/s，相差 %s（允许 %s，即理论值的 1%%）'
            % (_数(入口流量), _数(理论流量), _数(入口偏差), _数(入口限))),
        _检查稳定(行数, 数据['inletP'], 数据['outletT']),
        _检查残差(残差, 数据['outletT'][-残差窗口:], 入口温度),
    ]
    水侧检查, 水侧结果 = ([], {})
    if 信息.get('水侧'):
        水侧检查, 水侧结果 = _检查水侧(目录, 信息, 吸热量, 对照一维)
        检查 += [项 for 项 in 水侧检查 if 项['项目'] != '水出口温度与一维解一致']
    if 对照一维:
        压差偏差 = abs(多孔区压差 - 压降估算)
        压差限 = 一维容差 * 压降估算
        温度偏差 = abs(出口温度 - 温度估算)
        温度限 = 一维容差 * abs(入口温度 - 温度估算)
        检查 += [
            _项('压差与一维解一致', 压差偏差 < 压差限,
                '计算值 %s Pa，估算值 %s Pa，相差 %s Pa（允许偏差 %s Pa，即估算值的 3%%）'
                % (_数(多孔区压差), _数(压降估算), _数(压差偏差), _数(压差限))),
            _项('出口温度与一维解一致', 温度偏差 < 温度限,
                '计算值 %s K，估算值 %s K，相差 %s K（允许偏差 %s K，即估算温降 %s K 的 3%%）'
                % (_数(出口温度), _数(温度估算), _数(温度偏差), _数(温度限),
                   _数(abs(入口温度 - 温度估算)))),
        ]
        检查 += [项 for 项 in 水侧检查 if 项['项目'] == '水出口温度与一维解一致']

    结果 = {'通过': all(项['通过'] for 项 in 检查), '迭代步数': 行数, '进出口压差': 进出口压差,
            '多孔区压差': 多孔区压差,
            '入口流量': 入口流量, '出口流量': 出口流量, '出口温度': 出口温度, '吸热量': 吸热量,
            '多孔压降估算': 压降估算, '出口温度估算': 温度估算, **水侧结果, '检查': 检查}
    结果路径 = os.path.join(目录, 结果文件名)
    try:
        文本 = json.dumps(结果, ensure_ascii=False, allow_nan=False, indent=2) + '\n'
        with open(结果路径, 'w', encoding='utf-8', newline='\n') as f:
            f.write(文本)
    except (OSError, ValueError) as e:
        _清旧结果(目录)
        raise 结果错误('写不了 %s: %s' % (结果路径, e)) from e
    return dict(结果, 算例目录=目录, 结果路径=结果路径)


def _命令行(参数):
    '''返回退出码: 全部通过 0，有检查不通过或不能检查 1，用法不对 2。'''
    import sys
    选项 = [项 for 项 in 参数 if 项.startswith('--')]
    位置 = [项 for 项 in 参数 if not 项.startswith('--')]
    if len(位置) != 1 or len(选项) > 1 or any(项 != '--对照一维' for 项 in 选项):
        print(__doc__)
        return 2
    try:
        结果 = check_result(位置[0], 对照一维=bool(选项))
    except 结果错误 as e:
        print('不能检查: %s' % e, file=sys.stderr)
        return 1
    for 项 in 结果['检查']:
        print('%s  %s: %s' % ('通过' if 项['通过'] else '未通过', 项['项目'], 项['说明']))
    print('迭代 %d 步，进出口压差 %s Pa，入口流量 %s m³/s，出口温度 %s K，吸热量 %s W'
          % (结果['迭代步数'], _数(结果['进出口压差']), _数(结果['入口流量']),
             _数(结果['出口温度']), _数(结果['吸热量'])))
    print('结论: %s' % ('全部通过' if 结果['通过'] else '有检查项未通过'))
    print(结果['结果路径'])
    return 0 if 结果['通过'] else 1


if __name__ == '__main__':
    import sys
    sys.exit(_命令行(sys.argv[1:]))
