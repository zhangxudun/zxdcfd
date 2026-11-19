# -*- coding: utf-8 -*-
'''按 docs/05-写算例.md 写 OpenFOAM 算例: 模型.step + 第 2 步的边界和工况 + 第 4 步的 mesh.msh -> 算例目录。

用法:
    from generate_case import generate_case
    信息 = generate_case('模型.step', 边界, 网格, '算例目录', 工况=None, 步数=None)
    信息['算例目录']    # 里面有 0/、constant/、system/、Allrun、mesh.msh、算例信息.json

    边界: 第 2 步 generate_boundary 的返回值，或 标注.json 的路径。
    网格: 第 4 步 generate_mesh 的返回值，或 mesh.msh 的路径。
    工况: 边界是第 2 步的返回值时可以不给；边界是路径时必须给。
    步数: 正整数，不给时取 800。流向: 多孔区的主方向，默认 (0, 0, -1)。

命令行:
    python generate_case.py 模型.step 标注.json 工况.json mesh.msh 算例目录 [步数]

移植自 zhangxudun/cfd 的 openfoam算例.py，文件模板和算法不变:
    不可压，密度取入口温度下的值；常数有效粘度 = 湍流粘度 + 分子粘度（Sutherland 公式）；定常，simpleFoam。
    温度用 scalarTransport 在 simpleFoam 里一起解，多孔区的热源写在它的 fvOptions 里。
    多孔阻力 Δp/L = 粘性阻力系数·μ·v + 惯性阻力系数·½ρv²，μ 是分子粘度，横向乘以 横向倍数。
    多孔区恰好有两个交界面时给出一维估算（压降、出口温度），否则为空。
改动: 不再自己划网格和猜面类型；写之前核对边界名、分区名与 mesh.msh 一致；
    先写到临时目录，全部成功后才换入。任何一步失败都抛异常，并且不生成、不改动算例。

这一步只写文件，不运行 OpenFOAM。依赖: gmsh（读 STEP 量面积用）。
'''
import json
import math
import numbers
import os
import shutil
import uuid

from generate_boundary import 工况错误, 校验工况, 识别几何
from generate_mesh import _核对, _读边界, 多孔区名, 网格文件名


class 算例错误(ValueError):
    '''网格文件读不了、边界名或分区名与网格不一致、步数或流向不合法。'''


R气体 = 8.314462
默认步数 = 800
信息文件名 = '算例信息.json'
换行 = chr(10)
边界类型 = {'inlet': 'patch', 'outlet': 'patch', 'symmetry': 'symmetry', 'walls': 'wall'}
监测 = {'inletP': ('inlet', 'areaAverage', 'p', ''), 'inletFlow': ('inlet', 'sum', 'phi', ''),
        'outletFlow': ('outlet', 'sum', 'phi', ''),
        'outletT': ('outlet', 'weightedAverage', 'T', 'weightField     phi;')}
本步产物 = ('0', 'constant', 'system', 'Allrun', 网格文件名, 信息文件名)

_头 = '''FoamFile
{
    version     2.0;
    format      ascii;
    class       @类@;
    object      @名@;
}

'''
_块 = '''    @面@
    {
        @内容@
    }
'''
_温度方程 = '''    heatT
    {
        type            scalarTransport;
        libs            (solverFunctionObjects);
        writeControl    writeTime;
        field           T;
        D               @D@;
        nCorr           0;
        resetOnStartUp  no;
        bounded01       false;
        fvOptions
        {
@热源@        }
    }
'''
_热源 = '''            heat
            {
                type            scalarSemiImplicitSource;
                selectionMode   cellZone;
                cellZone        porous_zone;
                volumeMode      specific;
                sources
                {
@项@                }
            }
'''
_热源_温差 = '''                    T
                    {
                        explicit
                        {
                            type        exprField;
                            expression  "@式@";
                        }
                        implicit    @Sp@;
                    }
'''
_热源_定值 = '''                    T   (@Su@ 0);
'''
_监测项 = '''    @名@
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        writeFields     false;
        regionType      patch;
        name            @面@;
        operation       @操作@;
        @权@
        fields          (@场@);
    }
'''
_controlDict = '''application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         @步数@;
deltaT          1;
writeControl    timeStep;
writeInterval   @步数@;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;

functions
{
@温度@@监测@}
'''
_fvSchemes = '''ddtSchemes      { default steadyState; }
gradSchemes
{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss upwind;
    div(phi,T)      bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
'''
_fvSolution = '''solvers
{
    p
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.1;
    }
    T
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-10;
        relTol          0.001;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent      no;
    residualControl
    {
        p               1e-6;
        U               1e-6;
    }
}
relaxationFactors
{
    fields      { p 0.3; }
    equations   { U 0.7; T 1; }
}
'''
_多孔 = '''porosity
{
    type            explicitPorositySource;
    active          yes;

    explicitPorositySourceCoeffs
    {
        selectionMode   cellZone;
        cellZone        porous_zone;
        type            DarcyForchheimer;

        d   @d@;
        f   @f@;

        coordinateSystem
        {
            origin  (0 0 0);
            e1      @e1@;
            e2      @e2@;
        }
    }
}
'''
_场 = '''dimensions      @量纲@;
internalField   uniform @初值@;

boundaryField
{
@边界@}
'''
_改类型 = '''sed -i '/^ *@面@$/,/}/ s/type  *patch;/type @类@;/' constant/polyMesh/boundary
'''
_transportProperties = '''transportModel  Newtonian;
nu              [0 2 -1 0 0 0 0] @nu@;
'''
_turbulenceProperties = '''simulationType  laminar;
'''
_Allrun = '''#!/bin/bash
# 在装有 OpenFOAM 的环境(容器)里运行; 日志写在算例目录的 log.* 文件
cd "$(dirname "$0")"
trap 'chmod -R a+rwX .' EXIT
for f in /usr/lib/openfoam/openfoam*/etc/bashrc /opt/openfoam*/etc/bashrc; do
    if [ -f "$f" ]; then source "$f"; break; fi
done
set -eo pipefail
gmshToFoam mesh.msh > log.gmshToFoam 2>&1
@改类型@checkMesh > log.checkMesh 2>&1 || true
simpleFoam > log.simpleFoam 2>&1
'''


def _填(模板, **值):
    for 键, v in 值.items():
        模板 = 模板.replace('@%s@' % 键, str(v))
    return 模板


def _向量(v):
    return '(' + ' '.join('%.9g' % (x + 0.0) for x in v) + ')'


def _数(x):
    return '(%.9g)' % (x + 0.0)


def _写(路径, 文本, 可执行=False):
    os.makedirs(os.path.dirname(路径), exist_ok=True)
    with open(路径, 'w', encoding='utf-8', newline=换行) as f:
        f.write(文本)
    if 可执行:
        os.chmod(路径, 0o755)


def 粘度(T):
    '''Sutherland 公式, Pa·s'''
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def 流体参数(工况):
    物, T = 工况['物性'], 工况['入口']['温度']
    密度 = 物['工作压力'] * 物['摩尔质量'] / (R气体 * T)
    分子粘度 = 粘度(T)
    导热 = 物['湍流粘度'] * 物['比热'] / 物['湍流普朗特数'] + 分子粘度 * 物['比热'] / 0.7
    return {'密度': 密度, '分子粘度': 分子粘度, '有效运动粘度': (分子粘度 + 物['湍流粘度']) / 密度,
            '比热': 物['比热'], '热扩散系数': 导热 / (密度 * 物['比热'])}


def 垂直向量(v):
    '''返回 (单位化的 v, 与 v 垂直的单位向量)'''
    try:
        v = [float(x) for x in v]
    except (TypeError, ValueError):
        raise 算例错误('流向必须是三个数，实际为 %r' % (v,)) from None
    长 = math.sqrt(sum(x * x for x in v))
    if len(v) != 3 or not math.isfinite(长) or not 长 > 0:
        raise 算例错误('流向必须是三个有限数，且不能是零向量，实际为 %r' % (v,))
    v = [x / 长 for x in v]
    e = [0.0, 0.0, 0.0]
    e[min(range(3), key=lambda i: abs(v[i]))] = 1.0
    w = [v[1] * e[2] - v[2] * e[1], v[2] * e[0] - v[0] * e[2], v[0] * e[1] - v[1] * e[0]]
    长 = math.sqrt(sum(x * x for x in w))
    return v, [x / 长 for x in w]


def 阻力系数(工况, 流体):
    '''返回 (d, f), 第一个分量沿流向。OpenFOAM 的粘性项乘的是有效粘度, 所以 d 按 分子粘度/有效粘度 折算'''
    孔 = 工况['多孔']
    折算 = 流体['分子粘度'] / (流体['有效运动粘度'] * 流体['密度'])
    d, f, 倍 = 孔['粘性阻力系数'] * 折算, 孔['惯性阻力系数'], 孔['横向倍数']
    return (d, d * 倍, d * 倍), (f, f * 倍, f * 倍)


def 多孔压降估算(工况, 流体, 表观速度, 长度):
    '''均匀流过多孔区的压降, Pa'''
    孔 = 工况['多孔']
    return 长度 * (孔['粘性阻力系数'] * 流体['分子粘度'] * 表观速度
                   + 孔['惯性阻力系数'] * 0.5 * 流体['密度'] * 表观速度 ** 2)


def 出口温度估算(工况, 流体, 表观速度, 长度, 细分=2000):
    '''均匀流过多孔区后的温度, K'''
    换, T = 工况['多孔']['换热'], 工况['入口']['温度']
    G热容 = 流体['密度'] * 表观速度 * 流体['比热']
    ds = 长度 / 细分
    for n in range(细分):
        if 换['方式'] == '温差':
            Tw = 换['壁温顶部'] + (换['壁温底部'] - 换['壁温顶部']) * (n + 0.5) / 细分
            a = 换['h_av'] * ds / G热容
            T = (T * (1 - a / 2) + a * Tw) / (1 + a / 2)
        else:
            T = T + 换['体积热源'] * ds / G热容
    return T


def 热源项(工况, 流体, 流向, 多孔范围):
    '''温度方程(已除以 ρ·cp)在多孔区的源项; 多孔范围 = 多孔区沿流向的 (起点, 终点) 坐标'''
    换, 容 = 工况['多孔']['换热'], 流体['密度'] * 流体['比热']
    if 换['方式'] == '定热源':
        return _填(_热源_定值, Su='%.9g' % (换['体积热源'] / 容))
    n = 垂直向量(流向)[0]
    起, 终 = 多孔范围
    if not 终 > 起:
        raise 算例错误('多孔区沿流向的长度必须大于 0')
    s = '(%s*pos().x()+%s*pos().y()+%s*pos().z())' % (_数(n[0]), _数(n[1]), _数(n[2]))
    Tw = '(%s+%s*(%s-%s)/%s)' % (_数(换['壁温顶部']), _数(换['壁温底部'] - 换['壁温顶部']),
                                 s, _数(起), _数(终 - 起))
    k = 换['h_av'] / 容
    return _填(_热源_温差, 式='%s*%s' % (_数(k), Tw), Sp='%.9g' % -k)


def 写算例(算例目录, 工况, 面, 多孔范围, 流向=(0, 0, -1), 步数=默认步数):
    '''面: 网格里实际存在的边界名集合; 多孔范围: None 表示没有多孔区。返回流体参数。'''
    面 = set(面)
    流体 = 流体参数(工况)
    多孔, 热源 = '', ''
    if 多孔范围 is not None:
        e1, e2 = 垂直向量(流向)
        d, f = 阻力系数(工况, 流体)
        多孔 = _填(_多孔, d=_向量(d), f=_向量(f), e1=_向量(e1), e2=_向量(e2))
        热源 = _填(_热源, 项=热源项(工况, 流体, 流向, 多孔范围))

    def 写(相对路径, 类, 内容):
        路径 = os.path.join(算例目录, *相对路径.split('/'))
        _写(路径, _填(_头, 类=类, 名=os.path.basename(路径)) + 内容)

    def 边界(条件):
        return ''.join(_填(_块, 面=名, 内容=条件[名]) for 名 in sorted(面))

    写('system/controlDict', 'dictionary', _填(
        _controlDict, 步数=步数,
        温度=_填(_温度方程, D='%.9g' % 流体['热扩散系数'], 热源=热源),
        监测=''.join(_填(_监测项, 名=名, 面=面名, 操作=操作, 场=场, 权=权)
                     for 名, (面名, 操作, 场, 权) in 监测.items())))
    写('system/fvSchemes', 'dictionary', _fvSchemes)
    写('system/fvSolution', 'dictionary', _fvSolution)
    写('system/fvOptions', 'dictionary', 多孔)
    写('constant/transportProperties', 'dictionary',
       _填(_transportProperties, nu='%.9g' % 流体['有效运动粘度']))
    写('constant/turbulenceProperties', 'dictionary', _turbulenceProperties)
    T入 = '%.9g' % 工况['入口']['温度']
    写('0/U', 'volVectorField', _填(_场, 量纲='[0 1 -1 0 0 0 0]', 初值='(0 0 0)', 边界=边界({
        'inlet': 'type surfaceNormalFixedValue; refValue uniform %.9g;' % -工况['入口']['速度'],
        'outlet': 'type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);',
        'walls': 'type slip;' if 工况['壁面'] == '滑移' else 'type noSlip;',
        'symmetry': 'type symmetry;'})))
    写('0/p', 'volScalarField', _填(_场, 量纲='[0 2 -2 0 0 0 0]', 初值='0', 边界=边界({
        'inlet': 'type zeroGradient;',
        'outlet': 'type fixedValue; value uniform %.9g;' % (工况['出口']['表压'] / 流体['密度']),
        'walls': 'type zeroGradient;',
        'symmetry': 'type symmetry;'})))
    写('0/T', 'volScalarField', _填(_场, 量纲='[0 0 0 1 0 0 0]', 初值=T入, 边界=边界({
        'inlet': 'type fixedValue; value uniform %s;' % T入,
        'outlet': 'type inletOutlet; inletValue uniform %s; value uniform %s;' % (T入, T入),
        'walls': 'type zeroGradient;',
        'symmetry': 'type symmetry;'})))
    改 = ''.join(_填(_改类型, 面=名, 类=边界类型[名]) for 名 in ('symmetry', 'walls') if 名 in 面)
    _写(os.path.join(算例目录, 'Allrun'), _填(_Allrun, 改类型=改), 可执行=True)
    return 流体


def _校验步数(值):
    if 值 is None:
        return 默认步数
    if isinstance(值, bool) or not isinstance(值, numbers.Integral) or 值 <= 0:
        raise 算例错误('步数必须是正整数，实际为 %r' % (值,))
    return int(值)


def _读网格分组(路径):
    '''读 mesh.msh 开头的 $PhysicalNames，返回 (面分组名集合, 体分组名集合)。'''
    if not os.path.isfile(路径):
        raise 算例错误('找不到网格文件: %s' % 路径)
    if os.path.getsize(路径) == 0:
        raise 算例错误('网格文件是空的: %s' % 路径)
    组 = {2: set(), 3: set()}
    try:
        with open(路径, encoding='utf-8') as f:
            if f.readline().strip() != '$MeshFormat':
                raise 算例错误('网格文件不是 MSH 格式: %s' % 路径)
            在 = False
            for 行 in f:
                行 = 行.strip()
                if 行 in ('$EndPhysicalNames', '$Entities', '$Nodes', '$Elements'):
                    break
                if 行 == '$PhysicalNames':
                    在 = True
                elif 在:
                    段 = 行.split(None, 2)
                    if len(段) == 3 and 段[0] in ('2', '3'):
                        组[int(段[0])].add(段[2].strip().strip('"'))
    except UnicodeDecodeError:
        raise 算例错误('网格文件不是文本格式的 MSH: %s' % 路径) from None
    if not 组[2]:
        raise 算例错误('网格文件里没有面分组: %s' % 路径)
    return 组[2], 组[3]


def _是旧结果(目录, 名):
    if 名 in ('postProcessing', '结果.json') or 名.startswith('log.'):
        return True
    try:
        float(名)
    except ValueError:
        return False
    return os.path.isdir(os.path.join(目录, 名))


def _换入(临时, 目录):
    '''删掉旧算例和旧的求解结果，把临时目录里的新算例搬进来。目录里别的文件不动。'''
    for 名 in os.listdir(目录):
        路径 = os.path.join(目录, 名)
        if 路径 != 临时 and (名 in 本步产物 or _是旧结果(目录, 名)):
            if os.path.isdir(路径):
                shutil.rmtree(路径)
            else:
                os.remove(路径)
    for 名 in os.listdir(临时):
        os.replace(os.path.join(临时, 名), os.path.join(目录, 名))


def generate_case(step路径, 边界, 网格, 输出目录, 工况=None, 步数=None, 流向=(0, 0, -1)):
    '''核对边界与网格一致后写出算例，返回 算例信息.json 的内容，另加 算例目录、算例信息路径 两项。

    不改动传入的字典、模型.step、标注.json 和第 4 步的 mesh.msh。
    输出目录里只替换本步的产物和旧的求解结果（log.*、postProcessing、时间步目录、结果.json）。
    '''
    步数 = _校验步数(步数)
    n = 垂直向量(流向)[0]
    面标注, 多孔体 = _读边界(边界)
    if 工况 is None:
        if not isinstance(边界, dict) or '工况' not in 边界:
            raise 工况错误('没有工况: 边界不是第 2 步的返回值时，必须另外传入工况字典')
        工况 = 边界['工况']
    工况 = 校验工况(工况, bool(多孔体))
    try:
        网格路径 = os.fspath(网格['网格路径'] if isinstance(网格, dict) else 网格)
    except (KeyError, TypeError):
        raise 算例错误('网格必须是第 4 步的返回值或 mesh.msh 的路径，实际为 %r' % (网格,)) from None
    网格面组, 网格体组 = _读网格分组(网格路径)

    信息 = 识别几何(step路径)
    内部面 = {f['编号'] for f in 信息['面'] if len(f['相邻体']) == 2}
    外表面 = {f['编号'] for f in 信息['面']} - 内部面
    _核对(面标注, 多孔体, [v['编号'] for v in 信息['体']], 外表面, 内部面)
    面名 = {面标注[t] for t in 外表面}
    if 面名 != 网格面组:
        raise 算例错误('边界名与网格不一致: 只在边界里有 %s，只在网格里有 %s'
                       % (sorted(面名 - 网格面组), sorted(网格面组 - 面名)))
    if bool(多孔体) != (多孔区名 in 网格体组):
        raise 算例错误('分区名与网格不一致: 边界里%s多孔体，网格里%s %s 分区'
                       % ('有' if 多孔体 else '没有', '有' if 多孔区名 in 网格体组 else '没有', 多孔区名))
    if not {'inlet', 'outlet'} <= 面名:
        raise 算例错误('没有 inlet/outlet 面，请回到第 2 步在标注里指定面编号')

    多孔集 = set(多孔体)
    入口面积 = math.fsum(f['面积'] for f in 信息['面']
                         if f['编号'] in 外表面 and 面标注[f['编号']] == 'inlet')
    多孔范围 = None
    if 多孔集:                                   # 用多孔区各表面的中心估计它沿流向的起止位置
        s = [sum(a * b for a, b in zip(f['中心'], n)) for f in 信息['面'] if 多孔集 & set(f['相邻体'])]
        多孔范围 = (min(s), max(s))
    流体 = 流体参数(工况)
    流量 = 工况['入口']['速度'] * 入口面积
    交界 = [f['面积'] for f in 信息['面']
            if f['编号'] in 内部面 and len(多孔集 & set(f['相邻体'])) == 1]
    流速, 长度, 压降, 出口温度 = None, None, None, None
    if len(交界) == 2:                           # 多孔区一进一出两个交界面时才能估算
        截面 = 0.5 * sum(交界)
        流速 = 流量 / 截面
        长度 = math.fsum(v['体积'] for v in 信息['体'] if v['编号'] in 多孔集) / 截面
        压降 = 多孔压降估算(工况, 流体, 流速, 长度)
        出口温度 = 出口温度估算(工况, 流体, 流速, 长度)
    摘要 = {'边界': {名: 边界类型[名] for 名 in sorted(面名)}, '多孔区': 多孔区名 if 多孔集 else None,
            '步数': 步数, '流向': n, '监测量': sorted(监测),
            '入口面积': 入口面积, '入口速度': 工况['入口']['速度'], '入口温度': 工况['入口']['温度'],
            '出口表压': 工况['出口']['表压'], '密度': 流体['密度'], '比热': 流体['比热'],
            '分子粘度': 流体['分子粘度'], '有效运动粘度': 流体['有效运动粘度'],
            '理论入口体积流量': 流量, '理论入口质量流量': 流量 * 流体['密度'],
            '多孔范围': list(多孔范围) if 多孔范围 else None, '多孔区流速': 流速, '多孔区长度': 长度,
            '多孔压降估算': 压降, '出口温度估算': 出口温度}

    输出目录 = os.fspath(输出目录) or '.'
    os.makedirs(输出目录, exist_ok=True)
    临时 = os.path.join(输出目录, '.case.' + uuid.uuid4().hex + '.tmp')
    try:
        写算例(临时, 工况, 面名, 多孔范围, 流向, 步数)
        shutil.copyfile(网格路径, os.path.join(临时, 网格文件名))
        _写(os.path.join(临时, 信息文件名), json.dumps(摘要, ensure_ascii=False, indent=2) + 换行)
        _换入(临时, 输出目录)
    finally:
        shutil.rmtree(临时, ignore_errors=True)
    return dict(摘要, 算例目录=输出目录, 算例信息路径=os.path.join(输出目录, 信息文件名))


if __name__ == '__main__':
    import sys
    if len(sys.argv) not in (6, 7):
        print(__doc__)
        sys.exit(2)
    try:
        with open(sys.argv[3], encoding='utf-8') as f:
            输入工况 = json.load(f)
        结果 = generate_case(sys.argv[1], sys.argv[2], sys.argv[4], sys.argv[5], 工况=输入工况,
                               步数=int(sys.argv[6]) if len(sys.argv) == 7 else None)
    except (OSError, ValueError) as e:
        print('写算例不通过: %s' % e, file=sys.stderr)
        sys.exit(1)
    print('边界 %s，多孔区 %s，步数 %d' % (' / '.join(结果['边界']), 结果['多孔区'], 结果['步数']))
    print('理论入口流量 %.6g m³/s，多孔压降估算 %s Pa，出口温度估算 %s K'
          % (结果['理论入口体积流量'], 结果['多孔压降估算'], 结果['出口温度估算']))
    print(结果['算例目录'])
