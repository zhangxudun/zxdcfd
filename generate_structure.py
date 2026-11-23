# -*- coding: utf-8 -*-
'''按 docs/01-01-生成三维结构.md 生成 STEP（S1 直段管束、S2 变截面管束、S3 转弯入口管束）。

用法:
    from generate_structure import generate_step, 默认参数, 结构标注
    路径 = generate_step('S1', 默认参数('S1'), '输出目录', 'test_duct-01')
    标注 = 结构标注('S3', 参数, 路径)    # 第 2 步不能自动识别的结构给出标注，其余返回 None

依赖: gmsh（pip install gmsh）。坐标单位 m，STEP 文件声明的长度单位也是 m。
'''
import math
import numbers
import os
import re
import shutil
import tempfile
import uuid

import gmsh


class 参数错误(ValueError):
    '''参数缺失、非法、违反约束，或结构号不在清单里。'''


class 文件名错误(ValueError):
    '''输出文件名不合规。'''


_文件名规则 = re.compile('[A-Za-z0-9_-]+')
_毫米声明 = re.compile('SI_UNIT[(] *[.]MILLI[.] *, *[.]METRE[.] *[)]')
_米声明 = re.compile('SI_UNIT[(] *[$] *, *[.]METRE[.] *[)]')
_坐标规则 = re.compile('CARTESIAN_POINT *[(][^,]*, *[(]([^)]*)[)]')


def _长方体(x0, x1, y0, y1, z0, z1):
    gmsh.model.occ.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)


def _矩形(x0, x1, y0, y1, z):
    return [(x0, y0, z), (x1, y0, z), (x1, y1, z), (x0, y1, z)]


def _过渡段(底, 顶):
    '''两个平行矩形之间直线过渡的实体。底、顶各为 4 个角点，绕向和起点对应。'''
    occ = gmsh.model.occ
    线, 环 = [], []
    for 角点 in (底, 顶):
        点 = [occ.addPoint(*c) for c in 角点]
        边 = [occ.addLine(点[i], 点[(i + 1) % 4]) for i in range(4)]
        线 += 边
        环.append(occ.addCurveLoop(边))
    occ.addThruSections(环, makeSolid=True, makeRuled=True)
    occ.remove([(1, t) for t in 线], recursive=True)


def _建_S1(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    _长方体(-a, a, 0, b, h, h + li)
    _长方体(-a, a, 0, b, 0, h)
    _长方体(-a, a, 0, b, -lo, 0)


def _建_S2(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    ai, bi = p['入口半宽x'], p['入口深度y']
    ao, bo = p['出口半宽x'], p['出口深度y']
    # 开口在 y 方向居中，x 方向贴着对称面 x=0
    _过渡段(_矩形(0, a, 0, b, h), _矩形(0, ai, (b - bi) / 2, (b + bi) / 2, h + li))
    _长方体(0, a, 0, b, 0, h)
    _过渡段(_矩形(0, ao, (b - bo) / 2, (b + bo) / 2, -lo), _矩形(0, a, 0, b, 0))


def _建_S3(p):
    a, b = p['半宽x'], p['深度y']
    h, li, lo = p['管束高度'], p['入口段长度'], p['出口段长度']
    lh, hd = p['入口烟道长度'], p['入口烟道高度']
    顶 = h + li
    if _等高(hd, li):
        # 入口烟道与转向室等高，L 形退化成长方体
        _长方体(-a, a + lh, 0, b, h, 顶)
    else:
        # 在 y=0 平面画 L 形截面，沿 y 拉伸
        occ = gmsh.model.occ
        角点 = [(-a, h), (a, h), (a, 顶 - hd), (a + lh, 顶 - hd), (a + lh, 顶), (-a, 顶)]
        点 = [occ.addPoint(x, 0, z) for x, z in 角点]
        边 = [occ.addLine(点[i], 点[(i + 1) % 6]) for i in range(6)]
        面 = occ.addPlaneSurface([occ.addCurveLoop(边)])
        occ.extrude([(2, 面)], 0, b, 0)
    _长方体(-a, a, 0, b, 0, h)
    _长方体(-a, a, 0, b, -lo, 0)


def _约束_S2(p):
    for 小, 大 in (('入口半宽x', '半宽x'), ('入口深度y', '深度y'),
                   ('出口半宽x', '半宽x'), ('出口深度y', '深度y')):
        if p[小] > p[大]:
            raise 参数错误('违反约束: %s(%s) 不能大于 %s(%s)' % (小, p[小], 大, p[大]))


def _等高(甲, 乙):
    '''两个长度在 _等高容差 以内视为相等，避免算出来的 6.000000000001 走 L 形分支、生成细长面。'''
    return abs(甲 - 乙) <= _等高容差 * max(abs(甲), abs(乙))


def _约束_S3(p):
    if p['入口烟道高度'] > p['入口段长度'] and not _等高(p['入口烟道高度'], p['入口段长度']):
        raise 参数错误('违反约束: 入口烟道高度(%s) 不能大于 入口段长度(%s)'
                       % (p['入口烟道高度'], p['入口段长度']))


def _标注_S3(p):
    '''S3 的入口在侧面: 顶面会被第 2 步认成 inlet，要标成 walls；侧面入口要标成 inlet。'''
    a, b = p['半宽x'], p['深度y']
    顶 = p['管束高度'] + p['入口段长度']
    lh, hd = p['入口烟道长度'], p['入口烟道高度']
    return [((lh / 2, b / 2, 顶), 'walls'), ((a + lh, b / 2, 顶 - hd / 2), 'inlet')]


_等高容差 = 1e-9
_面心容差 = 1e-3         # m，按面心找面时允许的距离

# 新增结构: 在这里加一项（默认参数、约束函数、构造函数，第 2 步认不对的再加 标注），其余代码不用改。
_公共参数 = {'半宽x': 6.0, '深度y': 8.0, '管束高度': 3.2,
             '入口段长度': 3.0, '出口段长度': 4.0}

结构清单 = {
    'S1': {'名称': '直段管束', '参数': dict(_公共参数), '约束': None, '构造': _建_S1},
    'S2': {'名称': '变截面管束',
           '参数': dict(_公共参数, 入口半宽x=3.0, 入口深度y=4.0,
                        出口半宽x=2.0, 出口深度y=3.0),
           '约束': _约束_S2, '构造': _建_S2},
    'S3': {'名称': '转弯入口管束',
           '参数': dict(_公共参数, 入口段长度=6.0, 入口烟道长度=5.0, 入口烟道高度=4.0),
           '约束': _约束_S3, '构造': _建_S3, '标注': _标注_S3},
}


def _取结构(结构号):
    if not isinstance(结构号, str) or 结构号 not in 结构清单:
        raise 参数错误('未知结构号 %r，可选: %s' % (结构号, ', '.join(结构清单)))
    return 结构清单[结构号]


def 默认参数(结构号):
    '''返回该结构默认参数的副本。'''
    return dict(_取结构(结构号)['参数'])


def _校验参数(结构, 参数):
    if not isinstance(参数, dict):
        raise 参数错误('参数必须是字典')
    for 名 in 参数:
        if 名 not in 结构['参数']:
            raise 参数错误('未知参数: %s' % (名,))
    结果 = {}
    for 名 in 结构['参数']:
        if 名 not in 参数:
            raise 参数错误('缺少必填参数: %s' % 名)
        值 = 参数[名]
        if isinstance(值, bool) or not isinstance(值, numbers.Real):
            raise 参数错误('参数 %s 必须是数字，实际为 %r' % (名, 值))
        值 = float(值)
        if not math.isfinite(值):
            raise 参数错误('参数 %s 必须是有限数，实际为 %r' % (名, 值))
        if 值 <= 0:
            raise 参数错误('参数 %s 必须大于 0，实际为 %r' % (名, 值))
        结果[名] = 值
    if 结构['约束']:
        结构['约束'](结果)
    return 结果


def 结构标注(结构号, 参数, step路径):
    '''给出第 2 步要用的标注 {'面': {面编号: 类型}}；第 2 步能自动识别的结构返回 None。

    结构卡片按面心声明要改的面，这里用第 2 步的 识别几何 读回 STEP，按面心换成面编号，
    所以编号与第 2 步一致。step路径 必须是用同一组参数生成的文件。
    '''
    结构 = _取结构(结构号)
    p = _校验参数(结构, 参数)
    if not 结构.get('标注'):
        return None
    from generate_boundary import 识别几何      # 只有要标注的结构才用到第 2 步
    面 = 识别几何(step路径)['面']
    结果 = {}
    for 中心, 类 in 结构['标注'](p):
        命中 = [f['编号'] for f in 面 if math.dist(f['中心'], 中心) <= _面心容差]
        if len(命中) != 1:
            raise 参数错误('%s 标注: 面心在 %s 的面应恰好一个，实际 %d 个；'
                           'step路径 是否用同一组参数生成' % (结构号, 中心, len(命中)))
        结果[命中[0]] = 类
    return {'面': 结果}


def _校验文件名(文件名):
    if not isinstance(文件名, str) or not _文件名规则.fullmatch(文件名):
        raise 文件名错误('文件名 %r 不合规: 只能用英文字母、数字、下划线、连字符，'
                         '不能为空，不带扩展名' % (文件名,))


def _写出_单位m(临时文件):
    '''写出 STEP，保证声明的长度单位是 m，且坐标数值没有被换算。

    gmsh 导出时声明什么单位、数值是否换算，取决于版本和 Geometry.OCCTargetUnit。
    这里逐个尝试: 用文件里的最大坐标与模型最大坐标之比判断数值有没有被换算，
    数值没变的再把 mm 的单位声明改成 m。
    '''
    模型尺寸 = max(abs(v) for v in gmsh.model.getBoundingBox(-1, -1))
    for 单位 in ('M', 'MM', ''):
        gmsh.option.setString('Geometry.OCCTargetUnit', 单位)
        gmsh.write(临时文件)
        with open(临时文件, encoding='latin-1', newline='') as f:
            文本 = f.read()
        点组 = [组.split(',') for 组 in _坐标规则.findall(文本)]
        坐标 = [
            abs(float(值))
            for 点 in 点组 if len(点) == 3
            for 值 in 点
        ]
        if not 坐标 or not 0.5 < max(坐标) / 模型尺寸 < 2:
            continue
        文本 = _毫米声明.sub('SI_UNIT($,.METRE.)', 文本)
        if _米声明.search(文本):
            with open(临时文件, 'w', encoding='latin-1', newline='') as f:
                f.write(文本)
            return
    raise RuntimeError('STEP 写出失败: 得不到长度单位为 m 且数值未换算的文件')


def generate_step(结构号, 参数, 输出目录, 文件名):
    '''按参数生成 STEP，返回文件路径 <输出目录>/<文件名>.step。

    校验失败或生成失败时抛异常，并且不生成、不改动目标文件。
    '''
    结构 = _取结构(结构号)
    p = _校验参数(结构, 参数)
    _校验文件名(文件名)
    输出目录 = os.fspath(输出目录) or '.'
    目标 = os.path.join(输出目录, 文件名 + '.step')
    半成品 = os.path.join(输出目录, '.' + 文件名 + '.' + uuid.uuid4().hex + '.tmp')

    已初始化 = bool(gmsh.isInitialized())
    if 已初始化:
        原模型 = gmsh.model.getCurrent()
        原单位 = gmsh.option.getString('Geometry.OCCTargetUnit')
        原终端 = gmsh.option.getNumber('General.Terminal')
    else:
        gmsh.initialize()
    模型名 = 'gen_' + uuid.uuid4().hex
    临时目录 = tempfile.mkdtemp(prefix='genstep_')
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add(模型名)
        结构['构造'](p)
        gmsh.model.occ.synchronize()
        临时文件 = os.path.join(临时目录, 'model.step')
        _写出_单位m(临时文件)
        if not os.path.isfile(临时文件) or os.path.getsize(临时文件) == 0:
            raise RuntimeError('STEP 写出失败: 文件为空')
        os.makedirs(输出目录, exist_ok=True)
        shutil.copyfile(临时文件, 半成品)
        os.replace(半成品, 目标)
        return 目标
    finally:
        if os.path.exists(半成品):
            os.remove(半成品)
        shutil.rmtree(临时目录, ignore_errors=True)
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


if __name__ == '__main__':
    import sys
    号 = sys.argv[1] if len(sys.argv) > 1 else 'S1'
    目录 = sys.argv[2] if len(sys.argv) > 2 else '.'
    名 = sys.argv[3] if len(sys.argv) > 3 else 'model_' + 号
    print(generate_step(号, 默认参数(号), 目录, 名))
