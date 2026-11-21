# -*- coding: utf-8 -*-
'''网页界面（app.py）用到的工具: 三维预览用的面三角网格、默认工况、工况表单与字典的互转、读监测和残差、打包结果。

不依赖 streamlit，可以单独测试。依赖 gmsh。
'''
import copy
import io
import json
import os
import zipfile

import gmsh

from generate_boundary import _打开, 面类型

类型颜色 = {'inlet': '#2b7bba', 'outlet': '#d6604d', 'symmetry': '#7fbf7b', 'walls': '#bdbdbd',
            'interface': '#f4a582'}


def 默认工况():
    '''与 tests/test_end_to_end.py 的 多孔工况() 相同的一组“合理但人为”的数值，只用来起步，不代表真实设备。'''
    return {'入口': {'速度': 5.0, '温度': 600.0},
            '出口': {'表压': 0.0},
            '物性': {'工作压力': 101325.0, '摩尔质量': 0.029, '比热': 1100.0,
                     '湍流粘度': 1e-3, '湍流普朗特数': 0.85},
            '壁面': '滑移',
            '多孔': {'粘性阻力系数': 1e4, '惯性阻力系数': 10.0, '横向倍数': 100.0,
                     '换热': {'方式': '温差', 'h_av': 500.0, '壁温顶部': 500.0, '壁温底部': 450.0}}}


def 合并工况(底, 上):
    '''把 上 递归覆盖到 底 的副本上，返回 (合并结果, 不认识的键列表)。用于读入用户上传的 工况.json。'''
    结果, 不认识 = copy.deepcopy(底), []

    def 合(a, b, 路径):
        for 键, 值 in b.items():
            名 = '%s.%s' % (路径, 键) if 路径 else 键
            if 键 not in a and not (路径 == '多孔.换热'):
                不认识.append(名)
                continue
            if isinstance(值, dict) and isinstance(a.get(键), dict):
                合(a[键], 值, 名)
            else:
                a[键] = copy.deepcopy(值)
    if not isinstance(上, dict):
        raise ValueError('工况文件的最外层必须是字典')
    合(结果, 上, '')
    return 结果, 不认识


def 整理工况(工况, 有多孔体):
    '''去掉没有多孔体时用不到的 多孔 项，换热方式只留对应的键。'''
    结果 = copy.deepcopy(工况)
    if not 有多孔体:
        结果.pop('多孔', None)
        return 结果
    换 = 结果['多孔']['换热']
    留 = ('方式', 'h_av', '壁温顶部', '壁温底部') if 换['方式'] == '温差' else ('方式', '体积热源')
    结果['多孔']['换热'] = {k: 换[k] for k in 留 if k in 换}
    return 结果


def 面三角网格(step路径, 份数=12):
    '''返回 {面编号: {'点': [(x, y, z)], '三角形': [(i, j, k)]}}，面编号与第 2 步一致。份数 越大越细。'''
    结果 = {}
    with _打开(step路径):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(-1, -1)
        尺寸 = max(x1 - x0, y1 - y0, z1 - z0) / 份数
        gmsh.option.setNumber('Mesh.MeshSizeMax', 尺寸)
        gmsh.option.setNumber('Mesh.MeshSizeMin', 尺寸 / 4)
        gmsh.model.mesh.generate(2)
        for _, t in gmsh.model.getEntities(2):
            节点号, 坐标, _ = gmsh.model.mesh.getNodes(2, t, includeBoundary=True)
            位置 = {int(n): i for i, n in enumerate(节点号)}
            点 = [tuple(坐标[3 * i:3 * i + 3]) for i in range(len(节点号))]
            类型列, _, 连接列 = gmsh.model.mesh.getElements(2, t)
            三角 = []
            for 类, 连接 in zip(类型列, 连接列):
                if 类 != 2:          # 2 是三节点三角形
                    continue
                for i in range(0, len(连接), 3):
                    三角.append(tuple(位置[int(n)] for n in 连接[i:i + 3]))
            结果[int(t)] = {'点': 点, '三角形': 三角}
    return 结果


def 预览图(网格, 几何, 面类, 多孔体):
    '''用 plotly 画三维预览: 每个面按类型上色，鼠标停在面上显示面编号、类型和面积。返回 plotly Figure。'''
    import plotly.graph_objects as go
    面信息 = {f['编号']: f for f in 几何['面']}
    图 = go.Figure()
    已有图例 = set()
    for t, 片 in sorted(网格.items()):
        if not 片['三角形']:
            continue
        类 = 面类.get(t, 'walls')
        x, y, z = zip(*片['点'])
        i, j, k = zip(*片['三角形'])
        信息 = 面信息.get(t, {})
        图.add_trace(go.Mesh3d(
            x=x, y=y, z=z, i=i, j=j, k=k, color=类型颜色.get(类, '#999999'),
            opacity=0.35 if 类 == 'interface' else (0.6 if 类 == 'walls' else 0.95),
            name=类, legendgroup=类, showlegend=类 not in 已有图例, flatshading=True,
            hovertemplate='面 %d<br>%s<br>面积 %.4g m²<extra></extra>' % (t, 类, 信息.get('面积', 0))))
        已有图例.add(类)
    中心 = [f['中心'] for f in 几何['面'] if 面类.get(f['编号']) != 'interface']
    编号 = [f['编号'] for f in 几何['面'] if 面类.get(f['编号']) != 'interface']
    if 中心:
        cx, cy, cz = zip(*中心)
        图.add_trace(go.Scatter3d(x=cx, y=cy, z=cz, mode='text', text=[str(n) for n in 编号],
                                  textfont={'size': 12, 'color': '#111111'}, showlegend=False, hoverinfo='skip'))
    for v in 几何['体']:
        if v['编号'] in 多孔体:
            图.add_trace(go.Scatter3d(x=[v['中心'][0]], y=[v['中心'][1]], z=[v['中心'][2]], mode='markers+text',
                                      text=['多孔体 %d' % v['编号']], marker={'size': 5, 'color': '#b2182b'},
                                      showlegend=False, hoverinfo='skip'))
    图.update_layout(scene={'aspectmode': 'data', 'xaxis_title': 'x (m)', 'yaxis_title': 'y (m)',
                            'zaxis_title': 'z (m)'},
                     margin={'l': 0, 'r': 0, 't': 0, 'b': 0}, height=520,
                     legend={'orientation': 'h', 'y': 1.02})
    return 图


def 读历史(算例目录):
    '''读监测量和残差的逐步历史，供画图。返回 {'监测': {名: [值]}, '残差': {方程: [值]}}；读不到的为空。'''
    from check_result import _读监测, _读残差
    历史 = {'监测': {}, '残差': {}}
    for 名 in ('inletP', 'outletT', 'inletFlow', 'outletFlow'):
        try:
            历史['监测'][名] = _读监测(算例目录, 名)
        except Exception:
            pass
    try:
        历史['残差'] = _读残差(算例目录)
    except Exception:
        pass
    return 历史


def 打包(目录, 排除=('mesh.msh',)):
    '''把目录打成 zip（字节串），跳过 排除 里的文件名和 processor* 目录。'''
    缓冲 = io.BytesIO()
    with zipfile.ZipFile(缓冲, 'w', zipfile.ZIP_DEFLATED) as 包:
        for 根, 子目录, 文件 in os.walk(目录):
            子目录[:] = [d for d in 子目录 if not d.startswith('processor')]
            for 名 in 文件:
                if 名 in 排除:
                    continue
                路径 = os.path.join(根, 名)
                包.write(路径, os.path.relpath(路径, 目录))
    return 缓冲.getvalue()


def 转json(对象):
    return json.dumps(对象, ensure_ascii=False, indent=2, default=str)


可选面类型 = 面类型
