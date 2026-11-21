# -*- coding: utf-8 -*-
'''镜像冒烟测试（CI 在 docker 镜像里跑，不是 pytest 用例）:
内置 S1 → 提交到队列 → 后台 worker 在容器里用本机 OpenFOAM 真算一次 → 任务应“完成”，结论“通过”。
用法: docker run --rm zxdcfd:latest python tests/smoke_docker.py
'''
import os
import sys
import tempfile

sys.path.insert(0, '/app')
os.environ['ZXDCFD_WORKDIR'] = tempfile.mkdtemp()
os.environ['ZXDCFD_RUN_MODE'] = '本机'
from generate_structure import generate_step, 默认参数, 结构标注  # noqa: E402
import job_queue  # noqa: E402
import ui_tools  # noqa: E402
import worker  # noqa: E402

参数 = 默认参数('S1')
step = generate_step('S1', 参数, tempfile.mkdtemp(), 'model')
标注 = 结构标注('S1', 参数, step)
设置 = {'网格检查': False, '运行方式': '本机', '单元上限': 150000, '时限': 1800}
号 = job_queue.提交('冒烟', step, ui_tools.默认工况(), 标注, 设置, 说明='镜像冒烟')
worker.轮询间隔 = 2.0
worker.常驻(一次=True)
任务 = job_queue.读(号)
print(job_queue.日志(号, 末尾行数=40))
总结 = job_queue.总结(号) or {}
print('状态:', 任务['状态'], '结论:', 总结.get('结论'), 总结.get('说明'))
assert 任务['状态'] == '完成', 任务
assert 总结.get('结论') == '通过', 总结.get('说明')
print('镜像冒烟通过')
