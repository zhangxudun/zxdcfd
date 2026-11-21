# -*- coding: utf-8 -*-
'''网页界面 app.py 的测试（streamlit.testing.v1.AppTest，不开浏览器）。

UI-01  启动后提示先给几何，不报错。
UI-02  选内置 S1 生成几何: 自动识别出 1 个 inlet、1 个 outlet，多孔体 1 个，显示工况和设置。
UI-03  点“开始计算”: 用替身代替 auto_run.auto_run，确认界面把面类型、工况和设置原样传给它，并显示结论和结果表。
UI-04  内置 S3 生成几何: 套用结构卡片的标注，入口改到侧面，仍是 1 个 inlet。
需要 pip install streamlit plotly gmsh。
'''
import os
import sys

import pytest

仓库 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, 仓库)
pytest.importorskip('streamlit')
from streamlit.testing.v1 import AppTest  # noqa: E402

import auto_run  # noqa: E402


@pytest.fixture
def 界面(tmp_path, monkeypatch):
    monkeypatch.setenv('ZXDCFD_WORKDIR', str(tmp_path / '运行目录'))
    return AppTest.from_file(os.path.join(仓库, 'app.py'), default_timeout=120)


def 文字(at, 类):
    return '\n'.join(e.value for e in getattr(at, 类))


def 生成(at, 结构号):
    at.run()
    at.radio(key='几何来源').set_value('用内置结构生成').run()
    at.selectbox(key='结构号').set_value(结构号).run()
    at.button(key='生成几何').click().run()
    assert not at.exception, at.exception
    assert not at.error, 文字(at, 'error')


def test_UI_01_启动(界面):
    界面.run()
    assert not 界面.exception
    assert '先上传 STEP' in 文字(界面, 'info')


def test_UI_02_内置S1(界面):
    生成(界面, 'S1')
    assert '内置结构 S1' in 文字(界面, 'success')
    面 = 文字(界面, 'markdown')
    assert 'inlet 1 个' in 面 and 'outlet 1 个' in 面
    assert len(界面.multiselect[0].value) == 1          # 多孔体
    assert '还没有' not in 文字(界面, 'error')


def test_UI_03_开始计算(界面, tmp_path, monkeypatch):
    收到 = {}

    def 假运行(step路径, 工况, 输出目录, **参数):
        收到.update(参数, step路径=step路径, 工况=工况)
        参数['记']('[第 2 步] 替身')
        算例 = os.path.join(输出目录, '算例')
        os.makedirs(算例, exist_ok=True)
        结果 = {'进出口压差': 22.9, '出口温度': 498.8, '吸热量': 1.0e7, '入口流量': 160.0, '多孔压降估算': None,
                '出口温度估算': None, '检查': [{'项目': '流量守恒', '通过': True, '说明': '替身'}]}
        记录 = {'编号': 1, '层': '基准网格', '网格尺寸': 0.5, '单元数': 20000, '步数': 800, '完成步数': 800,
                '结论': {'结论': '通过'}, '改法': None, '异常': None, '结果': 结果, '算例目录': 算例}
        return {'结论': '通过', '说明': '替身说明', '最终': 记录, '尝试': [记录], '网格无关性': [],
                '输出目录': 输出目录}

    monkeypatch.setattr(auto_run, 'auto_run', 假运行)
    生成(界面, 'S1')
    界面.number_input(key='步数').set_value(1000)
    界面.checkbox(key='网格检查').uncheck()
    界面.button(key='开始计算').click().run()
    assert not 界面.exception, 界面.exception
    assert 收到['步数'] == 1000 and 收到['网格检查'] is False and 收到['网格尺寸'] is None
    assert sorted(收到['标注']['面'].values()).count('inlet') == 1
    assert len(收到['标注']['多孔体']) == 1
    assert 收到['工况']['入口']['速度'] == 5.0 and '多孔' in 收到['工况']
    assert '替身说明' in 文字(界面, 'success')
    assert any('22.9' in m.value for m in 界面.metric)


def test_UI_04_内置S3(界面):
    生成(界面, 'S3')
    面 = 文字(界面, 'markdown')
    assert 'inlet 1 个' in 面 and 'outlet 1 个' in 面
