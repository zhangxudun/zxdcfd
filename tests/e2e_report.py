# -*- coding: utf-8 -*-
'''端到端报告: 把 tests/test_end_to_end.py 走第 1–8 步时的实际结果记成文字。

这不是测试文件，不做断言。test_end_to_end.py 每算完一步就把这一步的结果记下来，测完写到一个文本文件，
所以每种结构只算一次，报告里的数值就是测试判定用的那一次。断言失败或某一步报错时，已经算出的步骤照样留在报告里；
报错交给第 8 步，记下返工结论。

报告路径取环境变量 E2E_REPORT，没有设置时写到仓库根目录的 e2e-report.txt。
CI 会把它发布到 ci-report 分支的 e2e-report.txt（PR 发布到 ci-report-pr 分支）。
'''
import collections
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from generate_structure import 结构清单  # noqa: E402
from check_result import _读残差, 残差窗口  # noqa: E402  与第 7 步“残差收敛”用同一个读法和窗口
from judge_result import judge_result, 判定错误  # noqa: E402

默认路径 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'e2e-report.txt')


def _点(坐标):
    return '(%s)' % ', '.join('%g' % (round(v, 6) + 0.0) for v in 坐标)


def _值(v):
    return '无' if v is None else '%.6g' % v


def _文(v):
    return json.dumps(v, ensure_ascii=False)


class 端到端报告:
    '''按结构分段记录，写出时按 结构清单 的顺序排。'''

    def __init__(self):
        self._段 = collections.OrderedDict()

    def 记录(self, 结构号):
        '''返回一个函数，调用它就往这种结构的段里加一行。'''
        行列 = self._段.setdefault(结构号, [])
        return 行列.append

    def 文本(self):
        结构号列 = [k for k in 结构清单 if k in self._段] + [k for k in self._段 if k not in 结构清单]
        return '\n'.join('\n'.join(self._段[k]) + '\n' for k in 结构号列)

    def 写出(self, 路径=None):
        路径 = 路径 or os.environ.get('E2E_REPORT') or 默认路径
        with open(路径, 'w', encoding='utf-8') as f:
            f.write(self.文本())
        return 路径


# ---------------------------------------------------------------- 各步

def 记输入(记, 结构号, 参数, 工况):
    记('%s %s %s %s' % ('=' * 20, 结构号, 结构清单[结构号]['名称'], '=' * 20))
    记('[输入] 结构参数: %s' % _文(参数))
    记('[输入] 工况: %s' % _文(工况))


def 记第1步(记, 模型, 标注):
    记('[第 1 步] 生成 %s，%d 字节' % (os.path.basename(模型), os.path.getsize(模型)))
    记('[输入] 标注: %s' % ('无' if 标注 is None else _文(标注)))


def 记第2步(记, 结果):
    几何 = 结果['几何']
    记('[第 2 步] 体 %d 个:' % len(几何['体']))
    for 体 in sorted(几何['体'], key=lambda v: v['编号']):
        记('  体 %d  体积 %.6g m³  中心 %s%s'
          % (体['编号'], 体['体积'], _点(体['中心']), '  多孔体' if 体['编号'] in 结果['多孔体'] else ''))
    记('[第 2 步] 面 %d 个:' % len(几何['面']))
    for 面 in sorted(几何['面'], key=lambda f: f['编号']):
        记('  面 %2d  %-9s  面积 %.6g m²  中心 %s  相邻体 %s'
          % (面['编号'], 结果['面'][面['编号']], 面['面积'], _点(面['中心']), 面['相邻体']))
    计数 = collections.Counter(结果['面'].values())
    记('[第 2 步] 面类型计数: %s' % ', '.join('%s %d' % 项 for 项 in sorted(计数.items())))
    记('[第 2 步] 多孔体: %s' % 结果['多孔体'])
    记('[第 2 步] 校验后的工况: %s' % _文(结果['工况']))
    记('[第 2 步] %s 内容:' % os.path.basename(结果['标注路径']))
    with open(结果['标注路径'], encoding='utf-8') as f:
        记(f.read().rstrip())


def 记第3步(记, 校验):
    记('[第 3 步] inlet 面 %s，outlet 面 %s' % (校验['inlet'], 校验['outlet']))


def 记第4步(记, 网格):
    记('[第 4 步] 生成 %s: 四面体 %d 个，节点 %d 个，网格尺寸 %s m'
      % (os.path.basename(网格['网格路径']), 网格['单元数'], 网格['节点数'], _值(网格['网格尺寸'])))
    记('[第 4 步] 网格体积 %s m³，几何体积 %s m³，相对偏差 %.3g%%'
      % (_值(网格['网格体积']), _值(网格['几何体积']), 100 * 网格['相对偏差']))
    for 体号, 体 in sorted(网格['体'].items()):
        记('  体 %d  四面体 %d 个  网格体积 %s m³  几何体积 %s m³  相对偏差 %.3g%%'
          % (体号, 体['单元数'], _值(体['网格体积']), _值(体['几何体积']), 100 * 体['相对偏差']))
    for 名, 组 in 网格['边界'].items():
        记('  %-9s  面 %s  网格面积 %s m²' % (名, 组['面'], _值(组['面积'])))
    记('[第 4 步] 体分组: %s' % _文(网格['体分组']))


def 记第5步(记, 算例):
    记('[第 5 步] 步数 %s，多孔区 %s' % (算例['步数'], 算例['多孔区']))
    记('[第 5 步] 边界: %s' % _文(算例['边界']))
    记('[第 5 步] 入口面积 %s m²，理论入口体积流量 %s m³/s，密度 %s kg/m³'
      % (_值(算例['入口面积']), _值(算例['理论入口体积流量']), _值(算例['密度'])))
    记('[第 5 步] 多孔区流速 %s m/s，多孔区长度 %s m，多孔压降估算 %s Pa，出口温度估算 %s K'
      % (_值(算例['多孔区流速']), _值(算例['多孔区长度']), _值(算例['多孔压降估算']), _值(算例['出口温度估算'])))
    if 算例.get('压力截面'):
        记('[第 5 步] 多孔区前后压力截面: %s' % '，'.join('%s 过 %s' % (名, _点(点)) for 名, 点 in 算例['压力截面'].items()))


def 记残差(记, 算例目录):
    try:
        残差 = _读残差(算例目录)
    except Exception as 错:          # 读不到残差不影响测试本身，第 7 步会给出正式结论
        记('[第 6 步] 残差读不了: %s' % 错)
        return
    记('[第 6 步] 残差（每步第一次求解的初始残差；fvSolution 里提前结束的条件是 p、U 都低于 1e-06）:')
    for 名, 值列 in 残差.items():
        首, 末, 末段 = 值列[0], 值列[-1], 值列[-残差窗口:]
        下降 = '%.1f' % math.log10(首 / 末) if 首 > 0 and 末 > 0 else '算不出'
        记('  %-3s 共 %d 步  第 1 步 %.3g  最后一步 %.3g  最后 %d 步里最大 %.3g 最小 %.3g  全程最大 %.3g  下降 %s 个数量级'
          % (名, len(值列), 首, 末, len(末段), max(末段), min(末段), max(值列), 下降))


def 记第6步(记, 求解):
    记('[第 6 步] 完成 %s 步（设定 %s 步）%s，监测量 %s'
      % (求解['完成步数'], 求解['步数'], '，提前收敛' if 求解['提前收敛'] else '',
         ' / '.join(sorted(求解['监测数据']))))
    记('[第 6 步] 转出的网格边界: %s' % _文(求解['边界']))
    记残差(记, 求解['算例目录'])


def 记第7步(记, 结果, 对照一维):
    记('[第 7 步] %s，写出 %s' % ('对照一维' if 对照一维 else '不对照一维', os.path.basename(结果['结果路径'])))
    for 项 in 结果['检查']:
        记('  %s  %s: %s' % ('通过' if 项['通过'] else '未通过', 项['项目'], 项['说明']))
    记('[第 7 步] 迭代 %d 步，进出口压差 %s Pa，多孔区压差 %s Pa，入口流量 %s m³/s，出口流量 %s m³/s，'
      '出口温度 %s K，吸热量 %s W'
      % (结果['迭代步数'], _值(结果['进出口压差']), _值(结果.get('多孔区压差')), _值(结果['入口流量']),
         _值(结果['出口流量']),
         _值(结果['出口温度']), _值(结果['吸热量'])))


def 记结论(记, 结论):
    for 项 in 结论['原因']:
        记('  %s，回到第 %s 步: %s。%s'
          % (项['现象'], '、'.join(str(n) for n in 项['回到步骤']), 项['处理'], 项['说明']))
    if 结论['结论'] == '通过':
        记('[第 8 步] 结论: 通过（%s）' % ('含一维对照' if 结论['对照一维'] else '未对照一维'))
    else:
        记('[第 8 步] 结论: 返工，回到第 %s 步' % '、'.join(str(n) for n in 结论['回到步骤']))


def 记中断(记, 错):
    '''某一步报错: 流程里的报错交给第 8 步记下返工结论；不是流程里的报错（程序出错）只记下异常。'''
    try:
        结论 = judge_result(错)
    except 判定错误:
        记('[程序出错] %s: %s' % (type(错).__name__, 错))
        return
    记('[中断] %s: %s' % (type(错).__name__, 错))
    记结论(记, 结论)
