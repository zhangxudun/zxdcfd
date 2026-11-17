# -*- coding: utf-8 -*-
'''按 docs/02-拿到边界条件.md 从 STEP 得到边界条件: 面类型、多孔体、工况数值。

用法:
    from generate_boundary import generate_boundary
    结果 = generate_boundary('模型.step', 工况, '输出目录', 标注=None)
    结果['标注路径']    # <输出目录>/标注.json，可直接交给后面的网格和算例步骤

命令行:
    python generate_boundary.py 模型.step 工况.json 输出目录 [标注.json]

面类型先自动猜: 两个体共用的面是 interface，位于 x=0 且是模型最小 x 的平面是 symmetry，
最高的水平面是 inlet，最低的水平面是 outlet，其余是 walls。几何里恰好有三个体时，
中间那个被猜成多孔体。猜错的用标注覆盖:
    标注 = {'面': {面编号: 'inlet' / 'outlet' / 'symmetry' / 'walls'}, '多孔体': [体编号]}

工况字典:
    {'入口': {'速度': m/s, '温度': K}, '出口': {'表压': Pa},
     '物性': {'工作压力': Pa, '摩尔质量': kg/mol, '比热': J/(kg·K),
              '湍流粘度': Pa·s, '湍流普朗特数': -},
     '壁面': '滑移' / '无滑移',
     '多孔': {'粘性阻力系数': 1/m², '惯性阻力系数': 1/m, '横向倍数': -,
              '换热': {'方式': '温差', 'h_av': W/(m³·K), '壁温顶部': K, '壁温底部': K}
                   或 {'方式': '定热源', '体积热源': W/m³}}}
    '多孔' 只在有多孔体时必填。

依赖: gmsh（pip install gmsh）。坐标单位 m，约定烟气大致沿 -z 流动。
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
    '物性': {'工作压力': _正, '摩尔质量': _正, '比热': _正,
             '湍流粘度': _非负, '湍流普朗特数': _正},
}
_多孔项 = {'粘性阻力系数': _非负, '惯性阻力系数': _非负, '横向倍数': _正}
_换热项 = {
    '温差': {'h_av': _非负, '壁温顶部': _正, '壁温底部': _正},
    '定热源': {'体积热源': _任意},
}


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


def _识别():
    '''对当前 gmsh 模型自动猜面类型和多孔体。'''
    occ = gmsh.model.occ
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(-1, -1)
    对角线 = math.dist((x0, y0, z0), (x1, y1, z1))
    容差 = 1e-4 * 对角线
    体 = [{'编号': int(t), '体积': occ.getMass(3, t),
           '中心': tuple(occ.getCenterOfMass(3, t))}
          for _, t in gmsh.model.getEntities(3)]
    面 = []
    for _, t in gmsh.model.getEntities(2):
        b = gmsh.model.getBoundingBox(2, t)
        相邻体 = [int(v) for v in gmsh.model.getAdjacencies(2, t)[0]]
        水平 = b[5] - b[2] < 容差
        if len(相邻体) == 2:
            类型 = 'interface'
        elif b[3] - b[0] < 容差 and abs(b[0]) < 容差 and abs(b[0] - x0) < 容差:
            类型 = 'symmetry'
        elif 水平 and abs(b[2] - z1) < 容差:
            类型 = 'inlet'
        elif 水平 and abs(b[2] - z0) < 容差:
            类型 = 'outlet'
        else:
            类型 = 'walls'
        面.append({'编号': int(t), '类型': 类型, '面积': occ.getMass(2, t),
                   '中心': tuple(occ.getCenterOfMass(2, t)), '相邻体': 相邻体})
    按高度 = sorted(体, key=lambda v: v['中心'][2])
    多孔体 = [按高度[1]['编号']] if len(体) == 3 else []
    return {'体': 体, '面': 面, '多孔体': 多孔体, '对角线': 对角线}


def 识别几何(step路径):
    '''只做自动识别，不看标注。返回 {'体': [...], '面': [...], '多孔体': [体编号], '对角线': m}。'''
    with _打开(step路径):
        return _识别()


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


def 校验工况(工况, 有多孔体):
    '''检查工况字典，返回数值都转成 float 的副本。不改动传入的字典。'''
    if not isinstance(工况, dict):
        raise 工况错误('工况必须是字典')
    结果 = copy.deepcopy(工况)
    for 组名, 项目 in _必填项.items():
        _取(结果, 组名, '')
        _校验一组(结果, 组名, 项目)
    壁面 = _取(结果, '壁面', '')
    if not isinstance(壁面, str) or 壁面 not in 壁面类型:
        raise 工况错误('工况的 壁面 只能是: %s，实际为 %r' % (' / '.join(壁面类型), 壁面))
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
    with _打开(step路径):
        信息 = _识别()
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
