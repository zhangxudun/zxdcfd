# -*- coding: utf-8 -*-
'''按 docs/03-校验.md 校验第 2 步得到的边界: 有 inlet 面和 outlet 面才能进入第 4 步。

用法:
    from validate_boundary import validate_boundary
    统计 = validate_boundary(结果)                  # 第 2 步 generate_boundary 的返回值
    统计 = validate_boundary('输出目录/标注.json')   # 或者 标注.json 的路径
    统计['inlet'], 统计['outlet']                    # 入口、出口的面编号，从小到大

命令行:
    python validate_boundary.py 标注.json
    通过时打印入口、出口的面编号并返回 0，不通过时打印原因并返回 1。

只看有没有，不看对不对: 入口、出口猜错了但存在，这一步照样通过。入口、出口可以不止一个。
这一步没有产物: 不写文件，不改传入的字典。不依赖 gmsh。
'''
import json
import numbers
import os
import sys


class 校验错误(ValueError):
    '''没有 inlet 面或 outlet 面，或者输入读不了、格式不对。'''


必须有的面类型 = ('inlet', 'outlet')


def _读输入(输入):
    '''输入可以是字典（第 2 步的返回值或标注），或 标注.json 的路径。返回 {面编号: 类型}。'''
    if not isinstance(输入, dict):
        try:
            路径 = os.fspath(输入)
        except TypeError:
            raise 校验错误('输入必须是第 2 步的返回值或 标注.json 的路径，实际为 %r' % (输入,)) from None
        try:
            with open(路径, encoding='utf-8') as f:
                输入 = json.load(f)
        except OSError as e:
            raise 校验错误('标注文件读不了: %s (%s)' % (路径, e)) from e
        except ValueError as e:
            raise 校验错误('标注文件不是合法的 JSON: %s (%s)' % (路径, e)) from e
        if not isinstance(输入, dict):
            raise 校验错误('标注文件的最外层必须是字典: %s' % 路径)
    if '面' not in 输入:
        raise 校验错误('输入里没有 面 这个键')
    if not isinstance(输入['面'], dict):
        raise 校验错误('输入的 面 必须是字典 {面编号: 类型}')
    return 输入['面']


def _转编号(值):
    if isinstance(值, numbers.Integral) and not isinstance(值, bool):
        return int(值)
    if isinstance(值, str):
        try:
            return int(值)
        except ValueError:
            pass
    raise 校验错误('面编号必须是整数，实际为 %r' % (值,))


def validate_boundary(输入):
    '''检查有没有 inlet 面和 outlet 面。

    输入: 第 2 步 generate_boundary 的返回值，或 标注.json 的路径。
    返回: {'inlet': [面编号], 'outlet': [面编号]}，面编号从小到大。
    缺 inlet 或 outlet 时抛 校验错误，信息里写明缺的是哪一种。
    '''
    面 = _读输入(输入)
    统计 = {类: sorted(_转编号(t) for t, 类型 in 面.items() if 类型 == 类)
            for 类 in 必须有的面类型}
    缺 = [类 for 类 in 必须有的面类型 if not 统计[类]]
    if 缺:
        raise 校验错误('没有 %s 面，请回到第 2 步在标注里指定面编号' % ' 面，也没有 '.join(缺))
    return 统计


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    try:
        结果 = validate_boundary(sys.argv[1])
    except 校验错误 as e:
        print('校验不通过: %s' % e, file=sys.stderr)
        sys.exit(1)
    print(', '.join('%s %d 个 %s' % (类, len(编号), 编号) for 类, 编号 in 结果.items()))
