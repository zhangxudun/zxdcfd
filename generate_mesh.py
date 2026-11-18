# -*- coding: utf-8 -*-
'''按 docs/04-划分网格.md 划分四面体网格: 模型.step + 第 2 步的边界 -> mesh.msh 和统计。

用法:
    from generate_mesh import generate_mesh
    统计 = generate_mesh('模型.step', 边界, '输出目录', 网格尺寸=None)
    统计['网格路径']    # <输出目录>/mesh.msh，MSH 2.2 文本格式

    边界: 第 2 步 generate_boundary 的返回值，或 标注.json 的路径。
    网格尺寸: 单位 m，不给时取模型包围盒对角线的 1/30。

命令行:
    python generate_mesh.py 模型.step 标注.json 输出目录 [网格尺寸]

网格里的分组: 面 inlet / outlet / symmetry / walls，体 porous_zone / fluid_zone。
内部面（interface）不进任何分组，也不写进 mesh.msh。
通过条件: 体单元全是四面体；每个体和总的网格体积与几何体积的相对偏差都不超过 1%。
任何一步失败都抛异常，并且不生成、不改动 mesh.msh。

这一步不重新猜面类型，也不检查有没有 inlet/outlet 面，那是第 2、3 步的事。
依赖: gmsh（pip install gmsh）。坐标单位 m。
'''
import contextlib
import json
import math
import numbers
import os
import uuid

import gmsh

from generate_boundary import _打开, _转编号, 标注错误, 面类型


class 网格错误(ValueError):
    '''网格尺寸不合法、出现非四面体单元，或网格体积与几何体积不符。'''


网格文件名 = 'mesh.msh'
体积容差 = 0.01          # “体积一致”的判定线: 相对偏差不超过 1%
默认尺寸比例 = 1.0 / 30   # 不给网格尺寸时，取包围盒对角线的这个比例
多孔区名, 流体区名 = 'porous_zone', 'fluid_zone'
_四面体, _三角形 = 4, 2   # gmsh 的单元类型号


def _校验尺寸(值):
    if 值 is None:
        return None
    if isinstance(值, bool) or not isinstance(值, numbers.Real):
        raise 网格错误('网格尺寸必须是数字，实际为 %r' % (值,))
    值 = float(值)
    if not math.isfinite(值) or 值 <= 0:
        raise 网格错误('网格尺寸必须是大于 0 的有限数，实际为 %r' % (值,))
    return 值


def _读边界(边界):
    '''边界可以是字典（第 2 步的返回值或标注），或 标注.json 的路径。返回 ({面编号: 类型}, [多孔体编号])。'''
    if not isinstance(边界, dict):
        try:
            路径 = os.fspath(边界)
        except TypeError:
            raise 标注错误('边界必须是第 2 步的返回值或 标注.json 的路径，实际为 %r' % (边界,)) from None
        try:
            with open(路径, encoding='utf-8') as f:
                边界 = json.load(f)
        except OSError as e:
            raise 标注错误('标注文件读不了: %s (%s)' % (路径, e)) from e
        except ValueError as e:
            raise 标注错误('标注文件不是合法的 JSON: %s (%s)' % (路径, e)) from e
        if not isinstance(边界, dict):
            raise 标注错误('标注文件的最外层必须是字典: %s' % 路径)
    if '面' not in 边界:
        raise 标注错误('边界里没有 面 这个键')
    if not isinstance(边界['面'], dict):
        raise 标注错误('边界的 面 必须是字典 {面编号: 类型}')
    面 = {}
    for 键, 类 in 边界['面'].items():
        t = _转编号(键, '面编号')
        if t in 面:
            raise 标注错误('面 %d 在边界里出现了两次' % t)
        if not isinstance(类, str) or 类 not in 面类型 + ('interface',):
            raise 标注错误('面 %d 的类型 %r 不支持，只能是: %s' % (t, 类, ' / '.join(面类型)))
        面[t] = 类
    原体 = 边界.get('多孔体', [])
    if not isinstance(原体, (list, tuple)):
        raise 标注错误('边界的 多孔体 必须是体编号的列表')
    多孔体 = [_转编号(v, '体编号') for v in 原体]
    if len(set(多孔体)) != len(多孔体):
        raise 标注错误('边界的 多孔体 里有重复的体编号')
    return 面, sorted(多孔体)


def _核对(面标注, 多孔体, 全部体, 外表面, 内部面):
    '''边界要和当前几何对得上: 每个外表面都有类型，没有多余的编号。'''
    不存在 = sorted(set(面标注) - 外表面 - 内部面)
    if 不存在:
        raise 标注错误('边界里有不存在的面编号: %s，标注和几何对不上' % 不存在)
    缺 = sorted(外表面 - set(面标注))
    if 缺:
        raise 标注错误('外表面没有类型: 面 %s，请回到第 2 步重新生成标注' % 缺)
    for t in sorted(内部面 & set(面标注)):
        if 面标注[t] != 'interface':
            raise 标注错误('面 %d 是两个体之间的内部面，不能设为边界' % t)
    for t in sorted(外表面):
        if 面标注[t] == 'interface':
            raise 标注错误('面 %d 是外表面，类型不能是 interface' % t)
    体不存在 = sorted(set(多孔体) - set(全部体))
    if 体不存在:
        raise 标注错误('体 %s 不存在' % 体不存在)


@contextlib.contextmanager
def _临时选项(选项):
    原 = {名: gmsh.option.getNumber(名) for 名 in 选项}
    try:
        for 名, 值 in 选项.items():
            gmsh.option.setNumber(名, 值)
        yield
    finally:
        for 名, 值 in 原.items():
            gmsh.option.setNumber(名, 值)


def _列表(x):
    return x.tolist() if hasattr(x, 'tolist') else list(x)


def _节点表():
    号, 坐标, _ = gmsh.model.mesh.getNodes()
    号, 坐标 = _列表(号), _列表(坐标)
    return {int(n): (坐标[3 * i], 坐标[3 * i + 1], 坐标[3 * i + 2]) for i, n in enumerate(号)}


def _四面体体积(a, b, c, d):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    wx, wy, wz = d[0] - a[0], d[1] - a[1], d[2] - a[2]
    return abs(ux * (vy * wz - vz * wy) - uy * (vx * wz - vz * wx) + uz * (vx * wy - vy * wx)) / 6.0


def _三角形面积(a, b, c):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    return math.sqrt((uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2) / 2.0


def _量(维, t, 点):
    '''返回 (单元数, 体积或面积)。体里只允许四面体，面上只允许三角形。'''
    应为, 每个, 算, 名 = ((_四面体, 4, _四面体体积, '四面体') if 维 == 3
                         else (_三角形, 3, _三角形面积, '三角形'))
    类型, _, 节点 = gmsh.model.mesh.getElements(维, t)
    数, 量 = 0, []
    for 型, 号 in zip(_列表(类型), 节点):
        if int(型) != 应为:
            raise 网格错误('%s %d 里有不是%s的单元（gmsh 单元类型 %d）'
                           % ('体' if 维 == 3 else '面', t, 名, int(型)))
        号 = _列表(号)
        for i in range(0, len(号), 每个):
            量.append(算(*(点[int(n)] for n in 号[i:i + 每个])))
        数 += len(号) // 每个
    return 数, math.fsum(量)


def _偏差(网格值, 几何值):
    return abs(网格值 - 几何值) / 几何值 if 几何值 > 0 else float('inf')


def _写网格(输出目录):
    目标 = os.path.join(输出目录, 网格文件名)
    半成品 = os.path.join(输出目录, '.mesh.' + uuid.uuid4().hex + '.tmp.msh')
    os.makedirs(输出目录, exist_ok=True)
    try:
        gmsh.write(半成品)
        os.replace(半成品, 目标)
    finally:
        if os.path.exists(半成品):
            os.remove(半成品)
    return 目标


def generate_mesh(step路径, 边界, 输出目录, 网格尺寸=None):
    '''划分四面体网格，检查通过后写出 <输出目录>/mesh.msh。

    返回:
        {'网格路径': 路径, '网格尺寸': m, '节点数': n, '单元数': 四面体总数,
         '网格体积': m³, '几何体积': m³, '相对偏差': -,
         '体': {体编号: {'单元数', '网格体积', '几何体积', '相对偏差'}},
         '边界': {分组名: {'面': [面编号], '面积': m²}},
         '体分组': {'porous_zone': [体编号], 'fluid_zone': [体编号]},   # 空的分组不出现
         '多孔体': [体编号]}

    不改动传入的字典、模型.step 和 标注.json。
    '''
    尺寸 = _校验尺寸(网格尺寸)
    面标注, 多孔体 = _读边界(边界)
    输出目录 = os.fspath(输出目录) or '.'
    with _打开(step路径):
        occ = gmsh.model.occ
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(-1, -1)
        if 尺寸 is None:
            尺寸 = math.dist((x0, y0, z0), (x1, y1, z1)) * 默认尺寸比例
        全部体 = sorted(int(t) for _, t in gmsh.model.getEntities(3))
        外表面, 内部面 = set(), set()
        for _, t in gmsh.model.getEntities(2):
            相邻体 = gmsh.model.getAdjacencies(2, t)[0]
            (内部面 if len(相邻体) == 2 else 外表面).add(int(t))
        _核对(面标注, 多孔体, 全部体, 外表面, 内部面)

        面组 = {类: sorted(t for t in 外表面 if 面标注[t] == 类) for 类 in 面类型}
        面组 = {类: 标签 for 类, 标签 in 面组.items() if 标签}
        体组 = {多孔区名: 多孔体, 流体区名: [t for t in 全部体 if t not in 多孔体]}
        体组 = {名: 标签 for 名, 标签 in 体组.items() if 标签}
        for 维, 组 in ((2, 面组), (3, 体组)):
            for 名, 标签 in 组.items():
                gmsh.model.setPhysicalName(维, gmsh.model.addPhysicalGroup(维, 标签), 名)
        几何体积 = {t: occ.getMass(3, t) for t in 全部体}

        with _临时选项({'Mesh.MeshSizeMax': 尺寸, 'Mesh.ElementOrder': 1,
                        'Mesh.MshFileVersion': 2.2, 'Mesh.Binary': 0, 'Mesh.SaveAll': 0}):
            gmsh.model.mesh.generate(3)
            点 = _节点表()
            体 = {}
            for t in 全部体:
                数, 体积 = _量(3, t, 点)
                if 数 == 0:
                    raise 网格错误('体 %d 没有划出四面体' % t)
                体[t] = {'单元数': 数, '网格体积': 体积, '几何体积': 几何体积[t],
                         '相对偏差': _偏差(体积, 几何体积[t])}
            总网格 = math.fsum(v['网格体积'] for v in 体.values())
            总几何 = math.fsum(几何体积.values())
            总偏差 = _偏差(总网格, 总几何)
            for 名, 项 in [('体 %d' % t, v) for t, v in 体.items()] + [
                    ('整个模型', {'网格体积': 总网格, '几何体积': 总几何, '相对偏差': 总偏差})]:
                if not 项['相对偏差'] <= 体积容差:
                    raise 网格错误(
                        '体积不符: %s 的网格体积 %.6g m³，几何体积 %.6g m³，相对偏差 %.3g%%，'
                        '超过 %g%%。请减小网格尺寸（当前 %g m），或修几何'
                        % (名, 项['网格体积'], 项['几何体积'], 100 * 项['相对偏差'],
                           100 * 体积容差, 尺寸))
            边界统计 = {类: {'面': 标签, '面积': math.fsum(_量(2, t, 点)[1] for t in 标签)}
                        for 类, 标签 in 面组.items()}
            路径 = _写网格(输出目录)
    return {'网格路径': 路径, '网格尺寸': 尺寸, '节点数': len(点),
            '单元数': sum(v['单元数'] for v in 体.values()),
            '网格体积': 总网格, '几何体积': 总几何, '相对偏差': 总偏差,
            '体': 体, '边界': 边界统计, '体分组': 体组, '多孔体': 多孔体}


if __name__ == '__main__':
    import sys
    if len(sys.argv) not in (4, 5):
        print(__doc__)
        sys.exit(2)
    try:
        结果 = generate_mesh(sys.argv[1], sys.argv[2], sys.argv[3],
                               float(sys.argv[4]) if len(sys.argv) == 5 else None)
    except ValueError as e:
        print('划分网格不通过: %s' % e, file=sys.stderr)
        sys.exit(1)
    print('四面体 %d 个，节点 %d 个，网格尺寸 %g m' % (结果['单元数'], 结果['节点数'], 结果['网格尺寸']))
    print('网格体积 %.6g m³，几何体积 %.6g m³，相对偏差 %.3g%%'
          % (结果['网格体积'], 结果['几何体积'], 100 * 结果['相对偏差']))
    print(', '.join('%s %.6g m²' % (类, v['面积']) for 类, v in 结果['边界'].items()))
    print(结果['网格路径'])
