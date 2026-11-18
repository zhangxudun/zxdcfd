# -*- coding: utf-8 -*-
'''端到端示例: 用默认参数把 S1、S2 各跑一遍第 1–2 步，打印每一步的实际结果。

这不是测试，不做断言，只是把 tests/test_end_to_end.py 走的那条链路的中间结果摆出来。
工况用 test_end_to_end.py 里的 多孔工况()，不加标注。
CI 会把输出发布到 ci-report 分支的 e2e-example.txt。

用法（在仓库根目录）:
    python tests/e2e_example.py [输出目录]
输出目录默认是 e2e-example，每种结构一个子目录，里面是 model.step 和 标注.json。
'''
import collections
import json
import os
import sys

_这里 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _这里)
sys.path.insert(0, os.path.dirname(_这里))
from generate_structure import generate_step, 默认参数, 结构清单  # noqa: E402
from generate_boundary import generate_boundary  # noqa: E402
from test_end_to_end import 多孔工况  # noqa: E402


def _点(坐标):
    return '(%s)' % ', '.join('%g' % (round(v, 6) + 0.0) for v in 坐标)


def 跑一个(结构号, 输出目录):
    参数 = 默认参数(结构号)
    目录 = os.path.join(输出目录, 结构号)
    print('=' * 20, 结构号, 结构清单[结构号]['名称'], '=' * 20)
    print('[输入] 结构参数: %s' % json.dumps(参数, ensure_ascii=False))
    print('[输入] 工况: %s' % json.dumps(多孔工况(), ensure_ascii=False))
    print('[输入] 标注: 无')

    模型 = generate_step(结构号, 参数, 目录, 'model')
    print('[第 1 步] 生成 %s，%d 字节' % (模型, os.path.getsize(模型)))

    结果 = generate_boundary(模型, 多孔工况(), 目录)
    几何 = 结果['几何']
    print('[第 2 步] 体 %d 个:' % len(几何['体']))
    for 体 in sorted(几何['体'], key=lambda v: v['编号']):
        print('  体 %d  体积 %.6g m³  中心 %s%s'
              % (体['编号'], 体['体积'], _点(体['中心']),
                 '  多孔体' if 体['编号'] in 结果['多孔体'] else ''))
    print('[第 2 步] 面 %d 个:' % len(几何['面']))
    for 面 in sorted(几何['面'], key=lambda f: f['编号']):
        print('  面 %2d  %-9s  面积 %.6g m²  中心 %s  相邻体 %s'
              % (面['编号'], 结果['面'][面['编号']], 面['面积'],
                 _点(面['中心']), 面['相邻体']))
    计数 = collections.Counter(结果['面'].values())
    print('[第 2 步] 面类型计数: %s'
          % ', '.join('%s %d' % 项 for 项 in sorted(计数.items())))
    print('[第 2 步] 多孔体: %s' % 结果['多孔体'])
    print('[第 2 步] 校验后的工况: %s' % json.dumps(结果['工况'], ensure_ascii=False))
    print('[第 2 步] %s 内容:' % 结果['标注路径'])
    with open(结果['标注路径'], encoding='utf-8') as f:
        print(f.read().rstrip())
    print()


if __name__ == '__main__':
    输出目录 = sys.argv[1] if len(sys.argv) > 1 else 'e2e-example'
    for 结构号 in 结构清单:
        跑一个(结构号, 输出目录)
