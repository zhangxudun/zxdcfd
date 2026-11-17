# -*- coding: utf-8 -*-
'''按 docs/01-01-生成三维结构.md 生成 STEP（S1 直段管束、S2 变截面管束）。

用法:
    from generate_structure import generate_step, 默认参数
    路径 = generate_step('S1', 默认参数('S1'), '输出目录', 'test_duct-01')

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


def _约束_S2(p):
    for 小, 大 in (('入口半宽x', '半宽x'), ('入口深度y', '深度y'),
                   ('出口半宽x', '半宽x'), ('出口深度y', '深度y')):
        if p[小] > p[大]:
            raise 参数错误('违反约束: %s(%s) 不能大于 %s(%s)' % (小, p[小], 大, p[大]))


# 新增结构: 在这里加一项（默认参数、约束函数、构造函数），其余代码不用改。
_公共参数 = {'半宽x': 6.0, '深度y': 8.0, '管束高度': 3.2,
             '入口段长度': 3.0, '出口段长度': 4.0}

结构清单 = {
    'S1': {'名称': '直段管束', '参数': dict(_公共参数), '约束': None, '构造': _建_S1},
    'S2': {'名称': '变截面管束',
           '参数': dict(_公共参数, 入口半宽x=3.0, 入口深度y=4.0,
                        出口半宽x=2.0, 出口深度y=3.0),
           '约束': _约束_S2, '构造': _建_S2},
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


def _校验文件名(文件名):
    if not isinstance(文件名, str) or not _文件名规则.fullmatch(文件名):
        raise 文件名错误('文件名 %r 不合规: 只能用英文字母、数字、下划线、连字符，'
                         '不能为空，不带扩展名' % (文件名,))


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
        gmsh.option.setString('Geometry.OCCTargetUnit', 'M')
        gmsh.model.add(模型名)
        结构['构造'](p)
        gmsh.model.occ.synchronize()
        临时文件 = os.path.join(临时目录, 'model.step')
        gmsh.write(临时文件)
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
