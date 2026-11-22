# -*- coding: utf-8 -*-
'''按 docs/02-拿到边界条件.md 从 STEP 得到边界条件: 面类型、多孔体、工况数值。

用法:
    from generate_boundary import generate_boundary
    结果 = generate_boundary('模型.step', 工况, '输出目录', 标注=None)
    结果['标注路径']    # <输出目录>/标注.json，可直接交给后面的网格和算例步骤

命令行:
    python generate_boundary.py 模型.step 工况.json 输出目录 [标注.json]

面类型先自动猜（流向默认 -z，可在工况里用 '流向' 改）: 两个体共用的面是 interface，位于 x=0 且是模型
最小 x 的平面是 symmetry（只在流向与 x=0 平面平行时才这样猜），垂直于流向、位于最上游的平面是 inlet，位于最下游的平面是 outlet，其余是 walls。
几何里恰好有三个体时，沿流向排在中间的那个被猜成多孔体。猜错的用标注覆盖:
    标注 = {'面': {面编号: 'inlet' / 'outlet' / 'symmetry' / 'walls'}, '多孔体': [体编号]}

工况字典（适用于任意单相牛顿流体，气体、液体都可以；带 * 的是可选项）:
    {'入口': {'速度': m/s, '温度': K}, '出口': {'表压': Pa},
     '物性': {'比热': J/(kg·K),
              密度二选一:   '密度': kg/m³（常数，液体或近似不可压的气体）
                            或 '工作压力': Pa + '摩尔质量': kg/mol（理想气体，按入口温度取值）
              分子粘度二选一: '分子粘度': Pa·s（常数）
                            或不填，按 Sutherland 公式随入口温度取值，系数 *'萨瑟兰参考粘度' Pa·s、
                            *'萨瑟兰参考温度' K、*'萨瑟兰常数' K 不填时用空气的值
              导热二选一:   '导热系数': W/(m·K) 或 '分子普朗特数': -
                            （都不填时只有理想气体可以，默认 0.7；常数密度的流体必须填）
              *'湍流粘度': Pa·s（常数，默认 0，即纯层流）, *'湍流普朗特数': -（默认 0.9）},
     *'流向': [x, y, z]（主流方向，默认 [0, 0, -1]）,
     '壁面': '滑移' / '无滑移',
     *'壁面换热': {'方式': '绝热'} 或 {'方式': '定壁温', '壁温': K}
                 或 {'方式': '定热流', '热流密度': W/m²（流入流体为正）}（默认绝热）,
     '多孔': {'粘性阻力系数': 1/m², '惯性阻力系数': 1/m, '横向倍数': -,
              '换热': {'方式': '温差', 'h_av': W/(m³·K), 壁温: 见下}
                   或 {'方式': '定热源', '体积热源': W/m³}}}
    温差换热的壁温: '壁温': K（均匀），或 '壁温入口侧' + '壁温出口侧': K（沿流向线性变化）；
        旧写法 '壁温顶部' / '壁温底部' 仍然认，含义同 入口侧 / 出口侧。
    '多孔' 只在有多孔体时必填。

依赖: gmsh（pip install gmsh）。坐标单位 m。
这一步不检查有没有 inlet/outlet 面，那是第 3 步校验的事。
'''
import contextlib
import copy
import json
import math
import numbers
import os
import uuid

import gmsh


class 几何错误(ValueError):
    '''STEP 不存在、读不进来，或里面没有实体。'''


class 标注错误(ValueError):
    '''标注格式不对，或面编号、体编号、面类型不合法。'''


class 工况错误(ValueError):
    '''工况字典缺项或数值不合法。'''


面类型 = ('inlet', 'outlet', 'symmetry', 'walls')
壁面类型 = ('滑移', '无滑移')
标注文件名 = '标注.json'

# (下限, 是否允许等于下限); 下限为 None 表示只要求是有限数
_正, _非负, _任意 = (0.0, False), (0.0, True), (None, False)
_必填项 = {
    '入口': {'速度': _正, '温度': _正},
    '出口': {'表压': _任意},
    '物性': {'比热': _正},
}
# 物性里可选的数值项。成对出现的（工作压力与摩尔质量）在 _校验物性 里另查
_物性可选项 = {'密度': _正, '工作压力': _正, '摩尔质量': _正, '分子粘度': _正,
               '萨瑟兰参考粘度': _正, '萨瑟兰参考温度': _正, '萨瑟兰常数': _非负,
               '导热系数': _正, '分子普朗特数': _正, '湍流粘度': _非负, '湍流普朗特数': _正}
默认流向 = (0.0, 0.0, -1.0)
默认湍流普朗特数 = 0.9
默认气体普朗特数 = 0.7
壁面换热方式 = {'绝热': {}, '定壁温': {'壁温': _正}, '定热流': {'热流密度': _任意}}
_多孔项 = {'粘性阻力系数': _非负, '惯性阻力系数': _非负, '横向倍数': _正}
_换热项 = {
    '温差': {'h_av': _非负},
    '定热源': {'体积热源': _任意},
}
_壁温别名 = {'壁温顶部': '壁温入口侧', '壁温底部': '壁温出口侧'}


@contextlib.contextmanager
def _打开(step路径):
    '''读入 STEP，多个实体先做 fragment。退出时恢复 gmsh 原来的状态。'''
    路径 = os.fspath(step路径)
    if not os.path.isfile(路径):
        raise 几何错误('找不到 STEP 文件: %s' % 路径)
    if os.path.getsize(路径) == 0:
        raise 几何错误('STEP 文件是空的: %s' % 路径)

    已初始化 = bool(gmsh.isInitialized())
    if 已初始化:
        原模型 = gmsh.model.getCurrent()
        原单位 = gmsh.option.getString('Geometry.OCCTargetUnit')
        原终端 = gmsh.option.getNumber('General.Terminal')
    else:
        gmsh.initialize()
    模型名 = 'bc_' + uuid.uuid4().hex
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
        gmsh.model.add(模型名)
        occ = gmsh.model.occ
        try:
            形体 = occ.importShapes(路径)
        except Exception as e:
            raise 几何错误('STEP 读入失败: %s (%s)' % (路径, e)) from e
        体 = [dt for dt in 形体 if dt[0] == 3]
        if not 体:
            raise 几何错误('STEP 文件里没有实体: %s' % 路径)
        if len(体) > 1:
            occ.fragment(体[:1], 体[1:])
        occ.synchronize()
        yield
    finally:
        if 已初始化:
            if 模型名 in gmsh.model.list():
                gmsh.model.setCurrent(模型名)
                gmsh.model.remove()
            if 原模型 in gmsh.model.list():
                gmsh.model.setCurrent(原模型)
            gmsh.option.setString('Geometry.OCCTargetUnit', 原单位)
            gmsh.option.setNumber('General.Terminal', 原终端)
        else:
            gmsh.finalize()


def 读流向(工况):
    '''工况里的 '流向'（不填为 -z），返回单位向量 (x, y, z)。工况不是字典时也返回默认值，留给 校验工况 报错。'''
    值 = 工况.get('流向', 默认流向) if isinstance(工况, dict) else 默认流向
    try:
        v = [float(x) for x in 值]
    except (TypeError, ValueError):
        raise 工况错误('工况的 流向 必须是三个数，实际为 %r' % (值,)) from None
    长 = math.sqrt(sum(x * x for x in v)) if len(v) == 3 else float('nan')
    if len(v) != 3 or not all(math.isfinite(x) for x in v) or not 长 > 0:
        raise 工况错误('工况的 流向 必须是三个有限数，且不能是零向量，实际为 %r' % (值,))
    return tuple(x / 长 for x in v)


def _投影范围(维, t, n, 轴):
    '''实体沿 n 的坐标范围。n 平行于坐标轴时用包围盒（精确），否则在参数域上取样估计。'''
    b = gmsh.model.getBoundingBox(维, t)
    if 轴 is not None:
        i, 号 = 轴
        return (b[i], b[i + 3]) if 号 > 0 else (-b[i + 3], -b[i])
    点 = []
    for _, 顶点 in gmsh.model.getBoundary([(维, t)], combined=False, recursive=True):
        点 += list(gmsh.model.getValue(0, 顶点, []))
    if 维 == 2:
        try:
            低, 高 = gmsh.model.getParametrizationBounds(2, t)
            参数 = []
            for i in range(7):
                for j in range(7):
                    参数 += [低[0] + (高[0] - 低[0]) * i / 6, 低[1] + (高[1] - 低[1]) * j / 6]
            点 += list(gmsh.model.getValue(2, t, 参数))
        except Exception:
            pass
    s = [点[k] * n[0] + 点[k + 1] * n[1] + 点[k + 2] * n[2] for k in range(0, len(点), 3)]
    if not s:
        s = [x * n[0] + y * n[1] + z * n[2] for x in (b[0], b[3]) for y in (b[1], b[4]) for z in (b[2], b[5])]
    return min(s), max(s)


def _识别(流向=默认流向):
    '''对当前 gmsh 模型自动猜面类型和多孔体。流向: 单位向量，inlet 在最上游，outlet 在最下游。'''
    occ = gmsh.model.occ
    n = tuple(流向)
    轴 = next(((i, 1 if n[i] > 0 else -1) for i in range(3) if abs(abs(n[i]) - 1) < 1e-12), None)
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(-1, -1)
    对角线 = math.dist((x0, y0, z0), (x1, y1, z1))
    容差 = 1e-4 * 对角线
    体 = [{'编号': int(t), '体积': occ.getMass(3, t),
           '中心': tuple(occ.getCenterOfMass(3, t))}
          for _, t in gmsh.model.getEntities(3)]
    范围 = {int(t): _投影范围(2, t, n, 轴) for _, t in gmsh.model.getEntities(2)}
    上游 = min(r[0] for r in 范围.values())
    下游 = max(r[1] for r in 范围.values())
    面 = []
    for _, t in gmsh.model.getEntities(2):
        b = gmsh.model.getBoundingBox(2, t)
        相邻体 = [int(v) for v in gmsh.model.getAdjacencies(2, t)[0]]
        低, 高 = 范围[int(t)]
        垂直 = 高 - 低 < 容差                      # 整个面沿流向的坐标都一样: 垂直于流向的平面
        if len(相邻体) == 2:
            类型 = 'interface'
        elif (abs(n[0]) < 1e-9 and b[3] - b[0] < 容差 and abs(b[0]) < 容差
              and abs(b[0] - x0) < 容差):           # 对称面要与流向平行
            类型 = 'symmetry'
        elif 垂直 and abs(低 - 上游) < 容差:
            类型 = 'inlet'
        elif 垂直 and abs(高 - 下游) < 容差:
            类型 = 'outlet'
        else:
            类型 = 'walls'
        面.append({'编号': int(t), '类型': 类型, '面积': occ.getMass(2, t),
                   '中心': tuple(occ.getCenterOfMass(2, t)), '相邻体': 相邻体})
    沿流向 = sorted(体, key=lambda v: sum(a * c for a, c in zip(v['中心'], n)))
    多孔体 = [沿流向[1]['编号']] if len(体) == 3 else []
    return {'体': 体, '面': 面, '多孔体': 多孔体, '对角线': 对角线, '流向': list(n)}


def 识别几何(step路径, 流向=默认流向):
    '''只做自动识别，不看标注。返回 {'体': [...], '面': [...], '多孔体': [体编号], '对角线': m, '流向': [...]}。'''
    with _打开(step路径):
        return _识别(流向)


def _转编号(值, 说明):
    if isinstance(值, bool):
        raise 标注错误('%s必须是整数，实际为 %r' % (说明, 值))
    if isinstance(值, numbers.Integral):
        return int(值)
    if isinstance(值, str) and 值.strip().lstrip('+-').isdigit():
        return int(值)
    raise 标注错误('%s必须是整数，实际为 %r' % (说明, 值))


def _读标注(标注):
    '''标注可以是 None、字典，或 JSON 文件路径。返回 (面标注, 多孔体或 None)。'''
    if 标注 is None:
        return {}, None
    if not isinstance(标注, dict):
        try:
            路径 = os.fspath(标注)
        except TypeError:
            raise 标注错误('标注必须是字典或 JSON 文件路径，实际为 %r' % (标注,)) from None
        try:
            with open(路径, encoding='utf-8') as f:
                标注 = json.load(f)
        except OSError as e:
            raise 标注错误('标注文件读不了: %s (%s)' % (路径, e)) from e
        except ValueError as e:
            raise 标注错误('标注文件不是合法的 JSON: %s (%s)' % (路径, e)) from e
        if not isinstance(标注, dict):
            raise 标注错误('标注文件的最外层必须是字典: %s' % 路径)
    for 键 in 标注:
        if 键 not in ('面', '多孔体'):
            raise 标注错误('标注里有未知的键: %r，只能有 面 和 多孔体' % (键,))

    原面 = 标注.get('面', {})
    if not isinstance(原面, dict):
        raise 标注错误('标注的 面 必须是字典 {面编号: 类型}')
    面标注 = {}
    for 键, 类 in 原面.items():
        t = _转编号(键, '面编号')
        if t in 面标注:
            raise 标注错误('面 %d 在标注里出现了两次' % t)
        if not isinstance(类, str) or 类 not in 面类型:
            raise 标注错误('面 %d 的类型 %r 不支持，只能是: %s' % (t, 类, ' / '.join(面类型)))
        面标注[t] = 类

    if '多孔体' not in 标注:
        return 面标注, None
    原体 = 标注['多孔体']
    if not isinstance(原体, (list, tuple)):
        raise 标注错误('标注的 多孔体 必须是体编号的列表')
    多孔体 = [_转编号(v, '体编号') for v in 原体]
    if len(set(多孔体)) != len(多孔体):
        raise 标注错误('标注的 多孔体 里有重复的体编号')
    return 面标注, 多孔体


def _套用标注(信息, 面标注, 多孔标注):
    面类 = {f['编号']: f['类型'] for f in 信息['面']}
    for t, 类 in 面标注.items():
        if t not in 面类:
            raise 标注错误('面 %d 不存在' % t)
        if 面类[t] == 'interface':
            raise 标注错误('面 %d 是两个体之间的内部面，不能设为边界' % t)
        面类[t] = 类
    if 多孔标注 is None:
        多孔体 = sorted(信息['多孔体'])
    else:
        全部体 = {v['编号'] for v in 信息['体']}
        不存在 = sorted(set(多孔标注) - 全部体)
        if 不存在:
            raise 标注错误('体 %s 不存在' % 不存在)
        多孔体 = sorted(多孔标注)
    return 面类, 多孔体


def _取(字典, 键, 路径):
    if not isinstance(字典, dict):
        raise 工况错误('工况的 %s 必须是字典' % 路径)
    if 键 not in 字典:
        raise 工况错误('工况缺少: %s.%s' % (路径, 键) if 路径 else '工况缺少: %s' % 键)
    return 字典[键]


def _数值(值, 名, 范围):
    下限, 可等于 = 范围
    if isinstance(值, bool) or not isinstance(值, numbers.Real):
        raise 工况错误('工况的 %s 必须是数字，实际为 %r' % (名, 值))
    值 = float(值)
    if not math.isfinite(值):
        raise 工况错误('工况的 %s 必须是有限数，实际为 %r' % (名, 值))
    if 下限 is not None and (值 < 下限 or (值 == 下限 and not 可等于)):
        raise 工况错误('工况的 %s 必须%s %g，实际为 %r'
                       % (名, '不小于' if 可等于 else '大于', 下限, 值))
    return 值


def _校验一组(结果, 组名, 项目):
    组 = 结果[组名]
    if not isinstance(组, dict):
        raise 工况错误('工况的 %s 必须是字典' % 组名)
    for 键, 范围 in 项目.items():
        名 = '%s.%s' % (组名, 键)
        组[键] = _数值(_取(组, 键, 组名), 名, 范围)


def _校验物性(物):
    '''物性里的可选项和二选一项。工作压力、摩尔质量要么都给，要么都不给。'''
    for 键, 范围 in _物性可选项.items():
        if 键 in 物:
            物[键] = _数值(物[键], '物性.' + 键, 范围)
    理想气体 = [键 for 键 in ('工作压力', '摩尔质量') if 键 in 物]
    if '密度' in 物:
        if 理想气体:
            raise 工况错误('工况的 物性 里 密度 与 %s 只能二选一: 给 密度 是常数密度，'
                           '给 工作压力 和 摩尔质量 是理想气体' % '、'.join(理想气体))
    elif len(理想气体) != 2:
        raise 工况错误('工况缺少: 物性.密度（常数密度），或 物性.工作压力 和 物性.摩尔质量（理想气体）'
                       + ('，现在只有 %s' % 理想气体[0] if 理想气体 else ''))
    萨瑟兰 = [键 for 键 in ('萨瑟兰参考粘度', '萨瑟兰参考温度', '萨瑟兰常数') if 键 in 物]
    if '分子粘度' in 物 and 萨瑟兰:
        raise 工况错误('工况的 物性 里 分子粘度 与 %s 只能二选一' % '、'.join(萨瑟兰))
    if '导热系数' in 物 and '分子普朗特数' in 物:
        raise 工况错误('工况的 物性 里 导热系数 与 分子普朗特数 只能二选一')
    if '密度' in 物 and '导热系数' not in 物 and '分子普朗特数' not in 物:
        raise 工况错误('工况缺少: 物性.导热系数 或 物性.分子普朗特数（常数密度的流体没有默认值）')


def _校验壁面换热(结果):
    if '壁面换热' not in 结果:
        return
    换 = 结果['壁面换热']
    if not isinstance(换, dict):
        raise 工况错误('工况的 壁面换热 必须是字典')
    方式 = _取(换, '方式', '壁面换热')
    if not isinstance(方式, str) or 方式 not in 壁面换热方式:
        raise 工况错误('工况的 壁面换热.方式 只能是: %s，实际为 %r' % (' / '.join(壁面换热方式), 方式))
    for 键, 范围 in 壁面换热方式[方式].items():
        换[键] = _数值(_取(换, 键, '壁面换热'), '壁面换热.' + 键, 范围)


def _校验壁温(换):
    '''温差换热的壁温: 均匀的 壁温，或沿流向线性的 壁温入口侧 + 壁温出口侧（旧名 壁温顶部 / 壁温底部）。'''
    for 旧, 新 in _壁温别名.items():
        if 旧 in 换 and 新 in 换:
            raise 工况错误('工况的 多孔.换热 里 %s 与 %s 是同一项，只能写一个' % (旧, 新))
    名 = {新: (旧 if 旧 in 换 else 新) for 旧, 新 in _壁温别名.items()}
    分端 = [名[新] for 新 in ('壁温入口侧', '壁温出口侧') if 名[新] in 换]
    if '壁温' in 换:
        if 分端:
            raise 工况错误('工况的 多孔.换热 里 壁温 与 %s 只能二选一' % '、'.join(分端))
        换['壁温'] = _数值(换['壁温'], '多孔.换热.壁温', _正)
        return
    if not 分端:
        raise 工况错误('工况缺少: 多孔.换热.壁温（均匀），或 多孔.换热.壁温入口侧 和 多孔.换热.壁温出口侧')
    for 新 in ('壁温入口侧', '壁温出口侧'):
        换[名[新]] = _数值(_取(换, 名[新], '多孔.换热'), '多孔.换热.' + 名[新], _正)


def 壁温两端(换):
    '''温差换热的 (入口侧壁温, 出口侧壁温)，K。换 必须已经过 校验工况。'''
    if '壁温' in 换:
        return 换['壁温'], 换['壁温']
    return (换.get('壁温入口侧', 换.get('壁温顶部')), 换.get('壁温出口侧', 换.get('壁温底部')))


def 校验工况(工况, 有多孔体):
    '''检查工况字典，返回数值都转成 float 的副本。不改动传入的字典。'''
    if not isinstance(工况, dict):
        raise 工况错误('工况必须是字典')
    结果 = copy.deepcopy(工况)
    for 组名, 项目 in _必填项.items():
        _取(结果, 组名, '')
        _校验一组(结果, 组名, 项目)
    _校验物性(结果['物性'])
    if '流向' in 结果:
        结果['流向'] = [float(x) for x in 读流向(结果)]
    壁面 = _取(结果, '壁面', '')
    if not isinstance(壁面, str) or 壁面 not in 壁面类型:
        raise 工况错误('工况的 壁面 只能是: %s，实际为 %r' % (' / '.join(壁面类型), 壁面))
    _校验壁面换热(结果)
    if not 有多孔体:
        return 结果

    孔 = _取(结果, '多孔', '')
    if not isinstance(孔, dict):
        raise 工况错误('工况的 多孔 必须是字典')
    for 键, 范围 in _多孔项.items():
        孔[键] = _数值(_取(孔, 键, '多孔'), '多孔.' + 键, 范围)
    换 = _取(孔, '换热', '多孔')
    方式 = _取(换, '方式', '多孔.换热')
    if not isinstance(方式, str) or 方式 not in _换热项:
        raise 工况错误('工况的 多孔.换热.方式 只能是: %s，实际为 %r'
                       % (' / '.join(_换热项), 方式))
    for 键, 范围 in _换热项[方式].items():
        换[键] = _数值(_取(换, 键, '多孔.换热'), '多孔.换热.' + 键, 范围)
    if 方式 == '温差':
        _校验壁温(换)
    return 结果


def _写标注(输出目录, 内容):
    目标 = os.path.join(输出目录, 标注文件名)
    半成品 = os.path.join(输出目录, '.' + 标注文件名 + '.' + uuid.uuid4().hex + '.tmp')
    os.makedirs(输出目录, exist_ok=True)
    try:
        with open(半成品, 'w', encoding='utf-8') as f:
            json.dump(内容, f, ensure_ascii=False, indent=2)
            f.write('\n')
        os.replace(半成品, 目标)
    finally:
        if os.path.exists(半成品):
            os.remove(半成品)
    return 目标


def generate_boundary(step路径, 工况, 输出目录, 标注=None):
    '''识别面类型和多孔体，套用标注，校验工况，写出 <输出目录>/标注.json。

    返回:
        {'标注路径': 路径,
         '面': {面编号: 类型},        # 全部面，内部面的类型是 interface
         '多孔体': [体编号],
         '工况': 校验后的工况副本,
         '几何': 自动识别的原始信息（体、面的编号、中心、面积等）}

    写出的 标注.json 只含外表面，每个外表面都有类型，格式与输入的标注相同。
    任何一步失败都抛异常，并且不生成、不改动 标注.json。
    '''
    面标注, 多孔标注 = _读标注(标注)
    流向 = 读流向(工况)
    with _打开(step路径):
        信息 = _识别(流向)
    面类, 多孔体 = _套用标注(信息, 面标注, 多孔标注)
    工况 = 校验工况(工况, bool(多孔体))
    输出目录 = os.fspath(输出目录) or '.'
    外表面 = {str(t): 类 for t, 类 in sorted(面类.items()) if 类 != 'interface'}
    路径 = _写标注(输出目录, {'面': 外表面, '多孔体': 多孔体})
    return {'标注路径': 路径, '面': 面类, '多孔体': 多孔体, '工况': 工况, '几何': 信息}


if __name__ == '__main__':
    import sys
    if len(sys.argv) not in (4, 5):
        print(__doc__)
        sys.exit(2)
    with open(sys.argv[2], encoding='utf-8') as f:
        输入工况 = json.load(f)
    结果 = generate_boundary(sys.argv[1], 输入工况, sys.argv[3],
                             sys.argv[4] if len(sys.argv) == 5 else None)
    统计 = {}
    for 类 in 结果['面'].values():
        统计[类] = 统计.get(类, 0) + 1
    print(', '.join('%s %d 个' % (类, 数) for 类, 数 in sorted(统计.items())))
    print('多孔体: %s' % 结果['多孔体'])
    print(结果['标注路径'])
