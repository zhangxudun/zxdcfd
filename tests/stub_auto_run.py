# -*- coding: utf-8 -*-
'''测试用的 auto_run 替身: 不划网格、不求解，按输入造一个“通过”的总结，并像 auto_run 一样写出 自动运行.json。

worker.py 在环境变量 ZXDCFD_AUTO_RUN=stub_auto_run:假运行 时用它（tests/ 目录要在 PYTHONPATH 里）。
ZXDCFD_STUB_SLEEP: 每写一行日志后等几秒，用来测试取消；ZXDCFD_STUB_FAIL=1 时抛异常，用来测试 出错。
'''
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import auto_run  # noqa: E402


def 假运行(step路径, 工况, 输出目录, 标注=None, 记=None, **设置):
    os.makedirs(输出目录, exist_ok=True)
    等 = float(os.environ.get('ZXDCFD_STUB_SLEEP') or 0)
    行列 = []
    for 行 in ('[第 2 步] 替身: 收到 STEP %s' % os.path.basename(step路径),
               '[第 3 步] 替身: 面 %s，多孔体 %s' % (json.dumps(标注['面'], ensure_ascii=False), 标注['多孔体']),
               '[第 1 次] 替身: 设置 %s' % json.dumps(设置, ensure_ascii=False, sort_keys=True)):
        行列.append(行)
        记(行)
        time.sleep(等)
    if os.environ.get('ZXDCFD_STUB_FAIL') == '1':
        raise RuntimeError('替身按要求出错')
    算例 = os.path.join(输出目录, '第01次', '算例')
    os.makedirs(算例, exist_ok=True)
    with open(os.path.join(算例, 'Allrun'), 'w') as f:
        f.write('# 替身\n')
    结果 = {'进出口压差': 22.9, '出口温度': 498.8, '吸热量': 1.0e7, '入口流量': 160.0, '多孔压降估算': None,
            '出口温度估算': None, '检查': [{'项目': '流量守恒', '通过': True, '说明': '替身'}]}
    记录 = {'编号': 1, '层': '基准网格', '网格尺寸': 0.5, '单元数': 20000, '步数': 设置.get('步数') or 800,
            '完成步数': 800, '结论': {'结论': '通过'}, '改法': None, '异常': None, '结果': 结果, '算例目录': 算例}
    总结 = {'结论': '通过', '说明': '替身说明', '最终': 记录, '尝试': [记录], '网格无关性': [], '边界': None,
            '输出目录': os.path.abspath(输出目录), '设置': 设置, '工况': 工况}
    行列.append('[结论] 通过。替身说明')
    记(行列[-1])
    auto_run._写出(总结, 行列)
    return 总结
