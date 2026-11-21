# -*- coding: utf-8 -*-
'''auto_run.py 的测试: 自动返工规则和网格无关性检查。

AR-01 至 AR-09 用替身代替第 4–7 步（generate_mesh / generate_case / run_case / check_result），
第 1–3 步和第 8 步（judge_result）用真的，所以不需要 OpenFOAM，几秒钟跑完。替身按网格尺寸、步数、时限
造出结果，用来确认每种出错都按 docs/09-网页界面与自动运行.md 第 3 节的规则自动修改，或者停下要人工。
AR-10 用容器真的运行 OpenFOAM: S3 从第 2 步自动跑到网格无关性检查结束，报告写到 AUTO_RUN_REPORT 指定的文件。
'''
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import auto_run  # noqa: E402
from generate_boundary import 识别几何  # noqa: E402
from generate_mesh import 网格错误  # noqa: E402
from generate_structure import generate_step, 默认参数, 结构标注  # noqa: E402
from run_case import 求解错误  # noqa: E402
import ui_tools  # noqa: E402

基本检查 = ('流量守恒', '入口流量', '监测值稳定', '残差收敛')


@pytest.fixture(scope='module')
def s1(tmp_path_factory):
    return generate_step('S1', 默认参数('S1'), str(tmp_path_factory.mktemp('S1')), 'S1')


class 替身:
    '''第 4–7 步的替身。参数都是函数，按 (网格尺寸, 步数, 时限) 决定行为。'''

    def __init__(self, monkeypatch, 压差=None, 吸热=None, 不通过=None, 网格报错=None, 求解报错=None,
                 单元系数=6 * 979.2):
        self.压差 = 压差 or (lambda h: 20.0 * (1 + 0.01 * h * h))
        self.吸热 = 吸热 or (lambda h: 1e7 * (1 + 0.01 * h * h))
        self.不通过 = 不通过 or (lambda h, n: [])
        self.网格报错 = 网格报错 or (lambda h: None)
        self.求解报错 = 求解报错 or (lambda h, n, t: None)
        self.单元系数 = 单元系数
        self.调用 = []
        for 名 in ('generate_mesh', 'generate_case', 'run_case', 'check_result'):
            monkeypatch.setattr(auto_run, 名, getattr(self, 名))

    def generate_mesh(self, 模型, 边界, 目录, 尺寸):
        self.调用.append(('网格', 尺寸))
        信息 = self.网格报错(尺寸)
        if 信息:
            raise 网格错误(信息)
        return {'网格尺寸': 尺寸, '单元数': int(self.单元系数 / 尺寸 ** 3), '相对偏差': 0.001}

    def generate_case(self, 模型, 边界, 网格, 目录, 步数=None):
        os.makedirs(目录, exist_ok=True)
        return {'算例目录': 目录, '步数': 步数, '网格尺寸': 网格['网格尺寸']}

    def run_case(self, 算例, 时限, 运行方式, 镜像):
        self.调用.append(('求解', 算例['网格尺寸'], 算例['步数'], 时限))
        信息 = self.求解报错(算例['网格尺寸'], 算例['步数'], 时限)
        if 信息:
            raise 求解错误(信息)
        return dict(算例, 完成步数=算例['步数'], 提前收敛=False)

    def check_result(self, 求解, 对照一维=False):
        h, n = 求解['网格尺寸'], 求解['步数']
        坏 = self.不通过(h, n)
        检查 = [{'项目': 名, '通过': 名 not in 坏, '说明': '替身'} for 名 in 基本检查]
        return {'通过': not 坏, '迭代步数': n, '进出口压差': self.压差(h), '入口流量': 160.0,
                '出口流量': 160.0, '出口温度': 500.0, '吸热量': self.吸热(h), '多孔压降估算': None,
                '出口温度估算': None, '检查': 检查, '算例目录': 求解['算例目录']}


def 跑(s1, 目录, **参数):
    行 = []
    参数.setdefault('网格检查', False)
    总结 = auto_run.auto_run(s1, ui_tools.默认工况(), str(目录), 记=行.append, **参数)
    return 总结, 行


def test_AR_01_首次通过_网格无关(s1, tmp_path, monkeypatch):
    替 = 替身(monkeypatch)
    总结, 行 = 跑(s1, tmp_path, 网格检查=True)
    assert 总结['结论'] == '通过', 总结['说明']
    assert len(总结['尝试']) == 3 and len(总结['网格无关性']) == 3   # 至少三层网格
    粗, 中, 细 = 总结['网格无关性']
    assert 中['网格尺寸'] == pytest.approx(粗['网格尺寸'] * auto_run.细化系数)
    assert 细['网格尺寸'] == pytest.approx(中['网格尺寸'] * auto_run.细化系数)
    assert 中['无关'] and 细['无关'] and 细['相对变化']['进出口压差'] < 0.03
    assert 总结['最终']['编号'] == 3                        # 取最细一层
    assert 粗['网格尺寸'] == pytest.approx(识别几何(s1)['对角线'] / 30)
    for 名 in (auto_run.总结文件名, auto_run.报告文件名):
        assert os.path.isfile(tmp_path / 名)
    assert any(l.startswith('[结论] 通过') for l in 行)
    assert 替.调用[0][0] == '网格'


def test_AR_02_残差不收敛_步数加倍(s1, tmp_path, monkeypatch):
    替身(monkeypatch, 不通过=lambda h, n: ['残差收敛'] if n < 1600 else [])
    总结, _ = 跑(s1, tmp_path)
    assert 总结['结论'] == '通过'
    第一, 第二 = 总结['尝试']
    assert '步数 800 → 1600' in 第一['改法']
    assert 第二['步数'] == 1600 and 第二['网格尺寸'] == 第一['网格尺寸']


def test_AR_03_网格错误_加密(s1, tmp_path, monkeypatch):
    默认尺寸 = 识别几何(s1)['对角线'] / 30
    替身(monkeypatch, 网格报错=lambda h: '网格体积与几何体积相差 2%' if h > 默认尺寸 * 0.9 else None)
    总结, _ = 跑(s1, tmp_path)
    assert 总结['结论'] == '通过'
    第一, 第二 = 总结['尝试']
    assert 第一['异常类'] == '网格错误' and '网格尺寸' in 第一['改法']
    assert 第二['网格尺寸'] == pytest.approx(默认尺寸 * auto_run.加密系数)


def test_AR_04_找不到docker_不重试(s1, tmp_path, monkeypatch):
    替 = 替身(monkeypatch, 求解报错=lambda h, n, t: "找不到 docker: 容器方式要先装好 docker")
    总结, _ = 跑(s1, tmp_path)
    assert 总结['结论'] == '需要人工'
    assert len(总结['尝试']) == 1 and '运行环境' in 总结['说明']
    assert sum(1 for c in 替.调用 if c[0] == '求解') == 1


def test_AR_05_没有入口_要人工(s1, tmp_path, monkeypatch):
    替 = 替身(monkeypatch)
    几何 = 识别几何(s1)
    标注 = {'面': {str(f['编号']): 'walls' for f in 几何['面'] if f['类型'] == 'inlet'}}
    总结, _ = 跑(s1, tmp_path, 标注=标注)
    assert 总结['结论'] == '需要人工' and 总结['尝试'] == []
    assert 'inlet' in 总结['说明']
    assert 替.调用 == []


def test_AR_06_网格无关性到单元上限(s1, tmp_path, monkeypatch):
    替身(monkeypatch, 压差=lambda h: 20.0 * (1 + 2 * h * h))
    总结, _ = 跑(s1, tmp_path, 网格检查=True, 单元上限=100000)
    assert 总结['结论'] == '网格无关性未确认'
    assert '超过上限' in 总结['说明']
    assert len(总结['网格无关性']) >= 2 and not 总结['网格无关性'][-1]['无关']
    assert all(r['单元数'] <= 100000 for r in 总结['尝试'])
    assert 总结['最终'] is not None


def test_AR_11_两层在容差内_不算无关(s1, tmp_path, monkeypatch):
    '''只加密一层时看不出是在收敛还是碰巧接近，不能判为无关。'''
    替身(monkeypatch)
    总结, _ = 跑(s1, tmp_path, 网格检查=True, 最多加密=1)
    assert 总结['结论'] == '网格无关性未确认'
    assert len(总结['网格无关性']) == 2 and 总结['网格无关性'][-1]['无关']


def test_AR_12_三层外推(s1, tmp_path, monkeypatch):
    '''压差 = 20(1 + 0.01 h²): 观测阶数 2，外推值 20；吸热量不随网格变化时 GCI 为 0；振荡时不外推。'''
    替身(monkeypatch, 吸热=lambda h: 1e7, 压差=lambda h: 20.0 * (1 + 0.01 * h * h))
    总结, 行 = 跑(s1, tmp_path, 网格检查=True)
    外推 = 总结['网格无关性'][-1]['外推']
    assert 外推['进出口压差']['阶数'] == pytest.approx(2, rel=1e-6)
    assert 外推['进出口压差']['外推值'] == pytest.approx(20, rel=1e-9)
    assert 外推['吸热量']['GCI'] == 0
    assert any('三层外推' in l for l in 行)
    三层 = [{'x': 1.0}, {'x': 2.0}, {'x': 1.5}]
    assert auto_run._外推(三层, 'x') is None


def test_AR_07_超时_时限加倍(s1, tmp_path, monkeypatch):
    替身(monkeypatch, 求解报错=lambda h, n, t: '超时: 求解超过时限 %g 秒' % t if t < 3600 else None)
    总结, _ = 跑(s1, tmp_path)
    assert 总结['结论'] == '通过'
    assert [r['时限'] for r in 总结['尝试']] == [1800, 3600]


def test_AR_08_返工次数用完(s1, tmp_path, monkeypatch):
    替身(monkeypatch, 不通过=lambda h, n: ['监测值稳定'])
    总结, _ = 跑(s1, tmp_path, 最多返工=2)
    assert 总结['结论'] == '未通过'
    assert [r['步数'] for r in 总结['尝试']] == [800, 1600, 3200]


def test_AR_09_入口流量_加密一次后要人工(s1, tmp_path, monkeypatch):
    替身(monkeypatch, 不通过=lambda h, n: ['入口流量'])
    总结, _ = 跑(s1, tmp_path)
    assert 总结['结论'] == '需要人工' and len(总结['尝试']) == 2
    assert '入口' in 总结['说明']


# ---------------------------------------------------------------- 真的运行 OpenFOAM

def docker可用():
    return shutil.which('docker') is not None and subprocess.run(
        ['docker', 'info'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def test_AR_10_S3_自动运行_网格无关性(tmp_path):
    if not docker可用():
        if os.environ.get('CI'):
            pytest.fail('CI 里不允许跳过 OF 类测试: 不能运行 docker')
        pytest.skip('不能运行 docker，OF 类测试跳过（跳过不算通过）')
    参数 = 默认参数('S3')
    step = generate_step('S3', 参数, str(tmp_path / '几何'), 'S3')
    行 = []
    总结 = auto_run.auto_run(step, ui_tools.默认工况(), str(tmp_path / '自动运行'), 标注=结构标注('S3', 参数, step),
                             网格检查=True, 单元上限=150000, 记=行.append)
    路径 = os.environ.get('AUTO_RUN_REPORT')
    if 路径:
        with open(tmp_path / '自动运行' / auto_run.报告文件名, encoding='utf-8') as f:
            报告 = f.read()
        with open(路径, 'w', encoding='utf-8') as f:
            f.write('S3 转弯入口管束，auto_run(网格检查=True, 单元上限=150000)\n\n' + 报告)
    assert 总结['结论'] in ('通过', '网格无关性未确认'), 总结['说明']
    assert len(总结['网格无关性']) >= 3, '至少要三层网格（加密两层）做比较'
    最终 = 总结['最终']
    assert 最终['结论']['结论'] == '通过'
    assert all(项['通过'] for 项 in 最终['结果']['检查'])
    assert 最终['单元数'] <= 150000
