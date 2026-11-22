# -*- coding: utf-8 -*-
'''按 docs/05-写算例.md 写 OpenFOAM 算例: 模型.step + 第 2 步的边界和工况 + 第 4 步的 mesh.msh -> 算例目录。

用法:
    from generate_case import generate_case
    信息 = generate_case('模型.step', 边界, 网格, '算例目录', 工况=None, 步数=None)
    信息['算例目录']    # 里面有 0/、constant/、system/、Allrun、mesh.msh、算例信息.json

    边界: 第 2 步 generate_boundary 的返回值，或 标注.json 的路径。
    网格: 第 4 步 generate_mesh 的返回值，或 mesh.msh 的路径。
    工况: 边界是第 2 步的返回值时可以不给；边界是路径时必须给。
    步数: 正整数，不给时取 800。流向: 主流方向，不给时取工况里的 '流向'，工况也没有时为 (0, 0, -1)。

命令行:
    python generate_case.py 模型.step 标注.json 工况.json mesh.msh 算例目录 [步数]

移植自 zhangxudun/cfd 的 openfoam算例.py，适用于任意单相牛顿流体（气体、液体）:
    不可压，密度是常数: 直接给的 密度，或理想气体在入口温度和工作压力下的值；
    分子粘度直接给，或 Sutherland 公式在入口温度下的值；定常，simpleFoam。
    层流（默认）: 常数有效粘度 = 湍流粘度（默认 0）+ 分子粘度；
        热扩散系数 = (分子导热 + 湍流粘度·cp/湍流普朗特数) / (ρ·cp)，分子导热 = 导热系数，或 μ·cp/分子普朗特数。
    湍流模型（k-ε、realizable k-ε、k-ω SST，RANS）: 粘度 = 分子粘度 + 模型算出的 νt；
        温度方程的扩散系数逐点为 α分子 + νt/湍流普朗特数；无滑移壁面用壁面函数（nutkWallFunction 等），
        入口按湍流强度和混合长度给 k、ε/ω，结束时写出壁面 y+。
    外壁: 绝热（零梯度）、定壁温或定热流，由工况的 壁面换热 决定。
    温度用 scalarTransport 在 simpleFoam 里一起解，多孔区的热源写在它的 fvOptions 里。
    温度方程另加 Euler 伪时间项，时间步按 伪时间步 取（温差换热时 Δt ≤ 1/k），防止换热汇项被显式处理后发散。
    多孔阻力 Δp/L = 粘性阻力系数·μ·v + 惯性阻力系数·½ρv²，μ 是分子粘度，横向乘以 横向倍数。
    多孔区恰好有两个交界面时给出一维估算（压降、出口温度），否则为空；同时在多孔区前后各加一个截面压力监测
    （porousInP、porousOutP），第 7 步用它们的差对照一维压降，不受进出口动压变化和转弯损失影响。
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

from generate_boundary import (工况错误, 校验工况, 识别几何, 读流向, 壁温两端, 默认气体普朗特数, 默认湍流普朗特数,
                               体积换热系数, 管束换热系数, 湍流模型名, 默认混合长度比例)
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
# 多孔区前后截面的平均压力（m²/s²），只在多孔区一进一出两个交界面时才加，供第 7 步对照一维压降。
# 截面垂直于流向，放在交界面外侧 截面偏移比例 × 多孔区长度 处，只切到紧挨交界面的流体单元，取单元值。
截面监测 = {'porousInP': '起', 'porousOutP': '终'}
截面偏移比例 = 1e-4
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
@扩散@
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
# 水侧耦合: 管内工质沿流向一维流动，多孔区沿流向分成 N 段，每个伪时间步:
#   1) 用当前烟气温度 T 和节点水温 Tn 算每段烟气传给水的热量 Q_b = Σ h_av·V·(T − Tw) 和 H_b = Σ h_av·V；
#   2) 从水入口一端往下游推: C水·(Tn下游 − Tn上游) = Q_b − H_b·(段平均水温新 − 段平均水温旧)
#      （把烟气温度暂时看成不变的线性化，比直接用 Q_b 推更稳），收敛后就是 C水·ΔTw = Q_b；
#   3) 温度方程加源项 k·(Tw − T)，−k·T 直接放进对角线（真正隐式）。
# 单元的水温在所在段两个节点之间按沿流向的位置线性插值。每步把水出口温度和换热量写进 postProcessing/waterSide.dat。
_水侧源 = '''            waterSide
            {
                type            scalarCodedSource;
                name            waterSide;
                selectionMode   cellZone;
                cellZone        porous_zone;
                fields          (T);

                codeInclude
                #{
                    #include <fstream>
                    #include <iomanip>
                    #include <cmath>
                #};

                codeCorrect
                #{
                #};

                codeConstrain
                #{
                #};

                codeAddSup
                #{
                    const scalar k = @k@;
                    const scalar hav = @hav@;
                    const scalar Cw = @Cw@;
                    const scalar TwIn = @TwIn@;
                    const label N = @N@;
                    const bool counter = @逆流@;
                    const vector n@n@;
                    const scalar s0 = @起@;
                    const scalar L = @长@;

                    static scalarField Tn(N + 1, TwIn);
                    static label lastIndex = -1;
                    static bool firstWrite = true;

                    const scalarField& V = mesh_.V();
                    const vectorField& C = mesh_.C();
                    const scalarField& T = eqn.psi().primitiveField();
                    const labelList& cs = cells_;

                    labelList bin(cs.size());
                    scalarField frac(cs.size());
                    forAll(cs, i)
                    {
                        const scalar x = (((C[cs[i]] & n) - s0)/L)*N;
                        const label b = min(max(label(std::floor(x)), label(0)), N - 1);
                        bin[i] = b;
                        frac[i] = min(max(x - b, scalar(0)), scalar(1));
                    }

                    const label idx = mesh_.time().timeIndex();
                    if (idx != lastIndex)
                    {
                        lastIndex = idx;
                        scalarField Q(N, Zero), H(N, Zero);
                        forAll(cs, i)
                        {
                            const label b = bin[i];
                            const scalar Tw = Tn[b] + frac[i]*(Tn[b + 1] - Tn[b]);
                            Q[b] += hav*V[cs[i]]*(T[cs[i]] - Tw);
                            H[b] += hav*V[cs[i]];
                        }
                        reduce(Q, sumOp<scalarField>());
                        reduce(H, sumOp<scalarField>());
                        const scalar Qsum = sum(Q);

                        scalarField Tnew(N + 1, TwIn);
                        for (label j = 0; j < N; ++j)
                        {
                            const label b = counter ? N - 1 - j : j;
                            const label up = counter ? b + 1 : b;
                            const label dn = counter ? b : b + 1;
                            const scalar a = H[b]/(2*Cw);
                            const scalar oldAvg = 0.5*(Tn[b] + Tn[b + 1]);
                            Tnew[dn] = (Tnew[up]*(1 - a) + (Q[b] + H[b]*oldAvg)/Cw)/(1 + a);
                        }
                        Tn = Tnew;

                        if (Pstream::master())
                        {
                            const fileName dir = mesh_.time().globalPath()/"postProcessing";
                            mkDir(dir);
                            std::ofstream os
                            (
                                (dir/"waterSide.dat").c_str(),
                                firstWrite ? std::ios::out : std::ios::app
                            );
                            os << std::setprecision(10);
                            if (firstWrite)
                            {
                                os << "# Time  TwOut  Qwater  Qgas" << std::endl;
                                firstWrite = false;
                            }
                            const scalar TwOut = counter ? Tn[0] : Tn[N];
                            os << mesh_.time().value() << ' ' << TwOut << ' '
                               << Cw*(TwOut - TwIn) << ' ' << Qsum << std::endl;

                            std::ofstream pf((dir/"waterProfile.dat").c_str());
                            pf << std::setprecision(10) << "# s  Tw" << std::endl;
                            for (label j = 0; j <= N; ++j)
                            {
                                pf << s0 + L*j/N << ' ' << Tn[j] << std::endl;
                            }
                        }
                    }

                    scalarField& src = eqn.source();
                    scalarField& dg = eqn.diag();
                    forAll(cs, i)
                    {
                        const label b = bin[i];
                        const scalar Tw = Tn[b] + frac[i]*(Tn[b + 1] - Tn[b]);
                        src[cs[i]] -= k*Tw*V[cs[i]];
                        dg[cs[i]] -= k*V[cs[i]];
                    }
                #};
            }
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
_截面监测项 = '''    @名@
    {
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        writeFields     false;
        regionType      sampledSurface;
        name            @名@;
        sampledSurfaceDict
        {
            type        plane;
            planeType   pointAndNormal;
            pointAndNormalDict
            {
                point   @点@;
                normal  @法向@;
            }
            interpolate false;
        }
        operation       areaAverage;
        fields          (p);
    }
'''
_controlDict = '''application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         @结束时间@;
deltaT          @时间步@;
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
_fvSchemes = '''ddtSchemes      { default steadyState; ddt(T) Euler; }
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
@湍流对流@}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
@壁面距离@'''
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
_turbulenceProperties_RAS = '''simulationType  RAS;

RAS
{
    RASModel        @模型@;
    turbulence      on;
    printCoeffs     on;
}
'''
_yPlus = '''    yPlus
    {
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        log             false;
    }
'''
Cmu = 0.09
湍流场 = {'kEpsilon': ('k', 'epsilon'), 'realizableKE': ('k', 'epsilon'), 'kOmegaSST': ('k', 'omega')}
# 湍流 + 无滑移壁面时，在 gmsh 四面体网格的壁面上挤出棱柱边界层（snappyHexMesh 只做 addLayers）。
# 纯四面体的近壁单元让壁面函数的摩擦偏低约 40%（docs/12），加 3 层后与光滑管关联式差 5–13%。
边界层 = {'层数': 3, '末层厚度比': 0.4, '扩展比': 1.2}
_加层命令 = 'snappyHexMesh -overwrite > log.snappyHexMesh 2>&1\n'
_snappyHexMeshDict = '''castellatedMesh false;
snap            false;
addLayers       true;

geometry { }

castellatedMeshControls
{
    maxLocalCells 1000000; maxGlobalCells 2000000; minRefinementCells 0; nCellsBetweenLevels 1;
    maxLoadUnbalance 0.1; features (); refinementSurfaces { } refinementRegions { }
    resolveFeatureAngle 30; locationInMesh (0 0 0); allowFreeStandingZoneFaces true;
}

snapControls { nSmoothPatch 3; tolerance 2; nSolveIter 30; nRelaxIter 5; }

addLayersControls
{
    relativeSizes   true;
    layers          { walls { nSurfaceLayers @层数@; } }
    expansionRatio  @扩展比@;
    finalLayerThickness @末层@;
    minThickness    0.02;
    nGrow           0;
    featureAngle    180;
    slipFeatureAngle 30;
    nRelaxIter      5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals  3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter      50;
}

meshQualityControls
{
    maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80;
    minVol 1e-15; minTetQuality -1e30; minArea -1; minTwist 0.02; minDeterminant 0.001;
    minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
    relaxed { maxNonOrtho 75; }
}

mergeTolerance 1e-6;
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
@改类型@@加层@checkMesh > log.checkMesh 2>&1 || true
@求解@
'''
_求解_普通 = 'simpleFoam > log.simpleFoam 2>&1'
# 水侧耦合的 scalarCodedSource 要在运行时编译，OpenFOAM 拒绝以 root 编译加载（dynamicCode::checkSecurity）。
# 容器里是 root 时改用 nobody 运行 simpleFoam（环境变量保留，HOME 指到 /tmp），先放开算例目录的写权限。
_求解_编译 = '''run_as=()
if [ "$(id -u)" = 0 ]; then
    chmod -R a+rwX .
    run_as=(setpriv --reuid=65534 --regid=65534 --clear-groups env HOME=/tmp)
fi
"${run_as[@]}" simpleFoam > log.simpleFoam 2>&1'''


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


def 粘度(T, 参考粘度=1.716e-5, 参考温度=273.15, 常数=110.4):
    '''Sutherland 公式, Pa·s。默认系数是空气的；别的气体在工况里给 萨瑟兰参考粘度/参考温度/常数'''
    return 参考粘度 * (T / 参考温度) ** 1.5 * (参考温度 + 常数) / (T + 常数)


def 流体参数(工况):
    '''入口温度下的常数物性。工况须已经过 校验工况。'''
    物, T = 工况['物性'], 工况['入口']['温度']
    if '密度' in 物:
        密度, 密度模型 = 物['密度'], '常数'
    else:
        密度, 密度模型 = 物['工作压力'] * 物['摩尔质量'] / (R气体 * T), '理想气体'
    if '分子粘度' in 物:
        分子粘度, 粘度模型 = 物['分子粘度'], '常数'
    else:
        分子粘度 = 粘度(T, 物.get('萨瑟兰参考粘度', 1.716e-5), 物.get('萨瑟兰参考温度', 273.15),
                        物.get('萨瑟兰常数', 110.4))
        粘度模型 = 'Sutherland'
    cp = 物['比热']
    if '导热系数' in 物:
        分子导热 = 物['导热系数']
    else:
        分子导热 = 分子粘度 * cp / 物.get('分子普朗特数', 默认气体普朗特数)
    湍流粘度 = 物.get('湍流粘度', 0.0)
    导热 = 分子导热 + 湍流粘度 * cp / 物.get('湍流普朗特数', 默认湍流普朗特数)
    return {'密度': 密度, '分子粘度': 分子粘度, '有效运动粘度': (分子粘度 + 湍流粘度) / 密度,
            '比热': cp, '热扩散系数': 导热 / (密度 * cp), '有效导热系数': 导热,
            '分子导热系数': 分子导热, '密度模型': 密度模型, '粘度模型': 粘度模型}


def 湍流参数(工况, 流体, 入口面积, 入口周长):
    '''入口当量直径 4A/P、入口雷诺数，选了湍流模型时还有入口的 k、ε 或 ω（按湍流强度 I 和混合长度 l）:
    k = 1.5·(U·I)²，ε = Cμ^0.75·k^1.5/l，ω = √k/(Cμ^0.25·l)，与 OpenFOAM 入口边界条件的公式相同。'''
    U = 工况['入口']['速度']
    当量直径 = 4 * 入口面积 / 入口周长 if 入口周长 > 0 else math.sqrt(4 * 入口面积 / math.pi)
    雷诺数 = U * 当量直径 * 流体['密度'] / 流体['分子粘度']
    模型 = 湍流模型名(工况)
    值 = {'模型': 工况.get('湍流', {}).get('模型', '层流'), 'OpenFOAM模型': 模型 or 'laminar',
          '入口当量直径': 当量直径, '入口雷诺数': 雷诺数}
    if 模型 is None:
        值['流态提示'] = ('入口雷诺数 %.3g 大于 4000，流动多半是湍流，层流结果只能看趋势，建议选湍流模型' % 雷诺数
                          if 雷诺数 > 4000 else None)
        return 值
    湍 = 工况['湍流']
    I = 湍['湍流强度']
    l = 湍.get('混合长度', 默认混合长度比例 * 当量直径)
    k = 1.5 * (U * I) ** 2
    值.update(湍流强度=I, 混合长度=l, 混合长度来源='工况' if '混合长度' in 湍 else '0.07×入口当量直径',
              k=k, epsilon=Cmu ** 0.75 * k ** 1.5 / l, omega=math.sqrt(k) / (Cmu ** 0.25 * l),
              湍流普朗特数=工况['物性'].get('湍流普朗特数', 默认湍流普朗特数),
              流态提示=('入口雷诺数 %.3g 小于 2300，流动可能是层流，RANS 湍流模型在低雷诺数下不可靠' % 雷诺数
                        if 雷诺数 < 2300 else None))
    return 值


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


def 水侧估算(工况, 流体, 表观速度, 长度, 截面):
    '''均匀流过多孔区、水侧耦合时的 ε-NTU 解: 返回 (烟气出口温度, 水出口温度)，K。
    C烟 = ρ·Q·cp，C水 = 水质量流量·水比热，UA = h_av·多孔区体积，NTU = UA/C小，Cr = C小/C大。'''
    换, T入 = 工况['多孔']['换热'], 工况['入口']['温度']
    C烟 = 流体['密度'] * 表观速度 * 截面 * 流体['比热']
    C水 = 换['水质量流量'] * 换['水比热']
    UA = 体积换热系数(换) * 截面 * 长度
    C小, C大 = min(C烟, C水), max(C烟, C水)
    NTU, Cr = UA / C小, C小 / C大
    if 换['水流动'] == '顺流':
        ε = (1 - math.exp(-NTU * (1 + Cr))) / (1 + Cr)
    elif abs(1 - Cr) < 1e-9:
        ε = NTU / (1 + NTU)
    else:
        e = math.exp(-NTU * (1 - Cr))
        ε = (1 - e) / (1 - Cr * e)
    热量 = ε * C小 * (T入 - 换['水入口温度'])           # 烟气传给水为正
    return T入 - 热量 / C烟, 换['水入口温度'] + 热量 / C水


def 出口温度估算(工况, 流体, 表观速度, 长度, 细分=2000, 截面=1.0):
    '''均匀流过多孔区后的温度, K'''
    换, T = 工况['多孔']['换热'], 工况['入口']['温度']
    if 换['方式'] == '水侧耦合':
        return 水侧估算(工况, 流体, 表观速度, 长度, 截面)[0]
    G热容 = 流体['密度'] * 表观速度 * 流体['比热']
    ds = 长度 / 细分
    for n in range(细分):
        if 换['方式'] == '温差':
            T入侧, T出侧 = 壁温两端(换)
            Tw = T入侧 + (T出侧 - T入侧) * (n + 0.5) / 细分
            a = 换['h_av'] * ds / G热容
            T = (T * (1 - a / 2) + a * Tw) / (1 + a / 2)
        else:
            T = T + 换['体积热源'] * ds / G热容
    return T


慢流停留时间 = 20.0     # s；平均停留时间超过它时放大伪时间步，见 伪时间步


def 伪时间步(工况, 流体, 停留时间=None):
    '''controlDict 的 deltaT，也是温度方程 Euler 伪时间项的步长 Δt（s）。其余方程是稳态，不受影响。

    scalarSemiImplicitSource 的 implicit 系数是负数时，OpenFOAM（fvm::SuSp）按显式处理，
    温差换热的汇项 -k·T（k = h_av/(ρ·cp)）用的是上一步的 T。对流很弱的单元里每步误差乘 (1 − k·Δt)，
    k·Δt > 2 时正负交替地发散。所以取 Δt ≤ 1/k：k ≤ 1 时 Δt = 1，否则取不超过 1/k 的最大的 2 的负整数次幂，
    让各步的时间值是精确的二进制小数。不用松弛因子: 松弛会作用到整条对流链上，强对流时温度每步只收敛千分之几。
    收敛后时间项为零，稳态解不变。没有多孔区或定热源时 Δt = 1。

    慢流（液体、低速）: 停留时间 = 流体域体积 / 入口体积流量（s）。超过 慢流停留时间 时，Δt = 1 的温度场
    要几百步才能流过一遍，流场先收敛、simpleFoam 提前停算，温度还没稳定。这时把 Δt 放大到不超过
    停留时间/2 的最大的 2 的整数次幂（仍受上面 1/k 的限制）。没有汇项时 Euler 隐式对任何 Δt 都稳定，
    Δt 越大越接近每步直接解稳态温度方程。停留时间不超过
    慢流停留时间 时（原来的气体算例都是这样）Δt 不变。
    '''
    放大 = 1.0
    if 停留时间 is not None and 停留时间 > 慢流停留时间:
        放大 = 2.0 ** math.floor(math.log2(停留时间 / 2))
    孔 = 工况.get('多孔')
    if not 孔 or 孔['换热']['方式'] not in ('温差', '水侧耦合'):
        return 放大
    k = 体积换热系数(孔['换热']) / (流体['密度'] * 流体['比热'])
    if not k > 0:
        return 放大
    return min(放大, 2.0 ** -math.ceil(math.log2(k)))     # 不超过 1/k 的最大的 2 的整数次幂


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
    T入侧, T出侧 = 壁温两端(换)
    Tw = '(%s+%s*(%s-%s)/%s)' % (_数(T入侧), _数(T出侧 - T入侧), s, _数(起), _数(终 - 起))
    k = 换['h_av'] / 容
    return _填(_热源_温差, 式='%s*%s' % (_数(k), Tw), Sp='%.9g' % -k)


def 水侧源项(工况, 流体, 流向, 多孔范围):
    '''水侧耦合时温度方程 fvOptions 里的 scalarCodedSource（运行时由 OpenFOAM 编译）。'''
    换 = 工况['多孔']['换热']
    起, 终 = 多孔范围
    if not 终 > 起:
        raise 算例错误('多孔区沿流向的长度必须大于 0')
    hav = 体积换热系数(换)
    return _填(_水侧源, k='%.12g' % (hav / (流体['密度'] * 流体['比热'])), hav='%.12g' % hav,
              Cw='%.12g' % (换['水质量流量'] * 换['水比热']), TwIn='%.12g' % 换['水入口温度'],
              N=换['分段数'], 逆流='true' if 换['水流动'] == '逆流' else 'false',
              n='(%s)' % ', '.join('%.12g' % (x + 0.0) for x in 垂直向量(流向)[0]), 起='%.12g' % 起, 长='%.12g' % (终 - 起))


def 水侧摘要(换):
    '''算例信息.json 里的 水侧 一项（没有水侧耦合时为 None）。'''
    if 换.get('方式') != '水侧耦合':
        return None
    值 = {'水入口温度': 换['水入口温度'], '水质量流量': 换['水质量流量'], '水比热': 换['水比热'],
          '水热容流量': 换['水质量流量'] * 换['水比热'], '水流动': 换['水流动'], '分段数': 换['分段数'],
          'h_av': 体积换热系数(换), '水饱和温度': 换.get('水饱和温度'),
          '监测文件': 'postProcessing/waterSide.dat', '水温分布文件': 'postProcessing/waterProfile.dat'}
    if '管束' in 换:
        管 = 管束换热系数(换['管束'])
        值.update(热阻=管['热阻'], 总传热系数=管['总传热系数'])
    return 值


def 壁面温度条件(工况, 流体):
    '''0/T 里 walls 的边界条件。定热流: 法向梯度 = q / 有效导热系数（q 流入流体为正）。
    湍流模型时壁面上的扩散系数是 α分子 + νt壁面/Prt（νt壁面由壁面函数给出，逐面不同），
    所以用 exprMixed 逐面算梯度 = q/(ρ·cp) / (α分子 + νt/Prt)，保证每个面的热流都是 q。'''
    换 = 工况.get('壁面换热', {'方式': '绝热'})
    if 换['方式'] == '定壁温':
        return 'type fixedValue; value uniform %.9g;' % 换['壁温']
    if 换['方式'] == '定热流' and 湍流模型名(工况):
        容 = 流体['密度'] * 流体['比热']
        Prt = 工况['物性'].get('湍流普朗特数', 默认湍流普朗特数)
        return ('type exprMixed; value uniform %.9g; valueExpr "0"; fractionExpr "0"; '
                'gradientExpr "%.12g/(%.12g + %.12g*nut)";'
                % (工况['入口']['温度'], 换['热流密度'] / 容, 流体['分子导热系数'] / 容, 1 / Prt))
    if 换['方式'] == '定热流':
        return 'type fixedGradient; gradient uniform %.9g;' % (换['热流密度'] / 流体['有效导热系数'])
    return 'type zeroGradient;'


def 截面位置(流向, 多孔范围):
    '''多孔区前后两个压力截面经过的点 {监测名: (x, y, z)}，分别在交界面外侧偏移一点的流体里。'''
    n = 垂直向量(流向)[0]
    起, 终 = 多孔范围
    偏移 = 截面偏移比例 * (终 - 起)
    位置 = {'起': 起 - 偏移, '终': 终 + 偏移}
    return {名: tuple(位置[端] * c for c in n) for 名, 端 in 截面监测.items()}


def _湍流文件(湍流, 工况, 边界, T入):
    '''选了湍流模型时的 0/k、0/epsilon 或 0/omega、0/nut: {相对路径: (类, 内容)}。'''
    模型 = 湍流['OpenFOAM模型']
    无滑移 = 工况['壁面'] == '无滑移'
    k, l = '%.9g' % 湍流['k'], '%.9g' % 湍流['混合长度']
    文件 = {'0/k': ('[0 2 -2 0 0 0 0]', k, {
        'inlet': 'type turbulentIntensityKineticEnergyInlet; intensity %.9g; value uniform %s;'
                 % (湍流['湍流强度'], k),
        'walls': 'type kqRWallFunction; value uniform %s;' % k if 无滑移 else 'type zeroGradient;'})}
    if 湍流场[模型][1] == 'epsilon':
        e = '%.9g' % 湍流['epsilon']
        文件['0/epsilon'] = ('[0 2 -3 0 0 0 0]', e, {
            'inlet': 'type turbulentMixingLengthDissipationRateInlet; mixingLength %s; value uniform %s;' % (l, e),
            'walls': 'type epsilonWallFunction; value uniform %s;' % e if 无滑移 else 'type zeroGradient;'})
    else:
        w = '%.9g' % 湍流['omega']
        文件['0/omega'] = ('[0 0 -1 0 0 0 0]', w, {
            'inlet': 'type turbulentMixingLengthFrequencyInlet; mixingLength %s; value uniform %s;' % (l, w),
            'walls': 'type omegaWallFunction; value uniform %s;' % w if 无滑移 else 'type zeroGradient;'})
    文件['0/nut'] = ('[0 2 -1 0 0 0 0]', '0', {
        'inlet': 'type calculated; value uniform 0;',
        'walls': 'type nutkWallFunction; value uniform 0;' if 无滑移 else 'type calculated; value uniform 0;'})
    结果 = {}
    for 路径, (量纲, 初值, 条件) in 文件.items():
        条件.setdefault('outlet', 'type inletOutlet; inletValue uniform %s; value uniform %s;' % (初值, 初值)
                        if 路径 != '0/nut' else 'type calculated; value uniform 0;')
        条件.setdefault('symmetry', 'type symmetry;')
        结果[路径] = _填(_场, 量纲=量纲, 初值=初值, 边界=边界(条件))
    return 结果


def _湍流求解设置(文本):
    '''fvSolution: k、ε、ω 与 U 用同样的求解器、松弛因子和残差控制。'''
    for 旧, 新 in (('    U\n    {', '    "(U|k|epsilon|omega)"\n    {'),
                   ('        U               1e-6;\n', '        U               1e-6;\n        "(k|epsilon|omega)" 1e-6;\n'),
                   ('{ U 0.7; T 1; }', '{ "(U|k|epsilon|omega)" 0.7; T 1; }')):
        assert 旧 in 文本, 旧
        文本 = 文本.replace(旧, 新)
    return 文本


def 要加边界层(工况, 面):
    '''选了湍流模型、壁面无滑移、网格里有 walls 时加棱柱边界层。'''
    return bool(湍流模型名(工况)) and 工况['壁面'] == '无滑移' and 'walls' in set(面)


def 写算例(算例目录, 工况, 面, 多孔范围, 流向=(0, 0, -1), 步数=默认步数, 截面=False, 停留时间=None, 湍流=None):
    '''面: 网格里实际存在的边界名集合; 多孔范围: None 表示没有多孔区。
    截面: True 时加多孔区前后截面的压力监测（要求有多孔范围）。
    湍流: 湍流参数() 的返回值；不给时按入口面积 1 m²、周长 4 m 估算（只影响入口 k、ε/ω 的默认混合长度）。
    返回流体参数。'''
    面 = set(面)
    流体 = 流体参数(工况)
    if 湍流 is None:
        湍流 = 湍流参数(工况, 流体, 1.0, 4.0)
    模型 = 湍流模型名(工况)
    if 模型:
        扩散 = '        alphaD          %.12g;\n        alphaDt         %.12g;' % (
            流体['分子导热系数'] / (流体['密度'] * 流体['比热']) / (流体['分子粘度'] / 流体['密度']),
            1 / 湍流['湍流普朗特数'])
    else:
        扩散 = '        D               %.9g;' % 流体['热扩散系数']
    多孔, 热源 = '', ''
    if 多孔范围 is not None:
        e1, e2 = 垂直向量(流向)
        d, f = 阻力系数(工况, 流体)
        多孔 = _填(_多孔, d=_向量(d), f=_向量(f), e1=_向量(e1), e2=_向量(e2))
        if 工况['多孔']['换热']['方式'] == '水侧耦合':
            热源 = 水侧源项(工况, 流体, 流向, 多孔范围)
        else:
            热源 = _填(_热源, 项=热源项(工况, 流体, 流向, 多孔范围))

    def 写(相对路径, 类, 内容):
        路径 = os.path.join(算例目录, *相对路径.split('/'))
        _写(路径, _填(_头, 类=类, 名=os.path.basename(路径)) + 内容)

    def 边界(条件):
        return ''.join(_填(_块, 面=名, 内容=条件[名]) for 名 in sorted(面))

    Δt = 伪时间步(工况, 流体, 停留时间)
    写('system/controlDict', 'dictionary', _填(
        _controlDict, 步数=步数, 时间步='%.9g' % Δt, 结束时间='%.9g' % (步数 * Δt),
        温度=_填(_温度方程, 扩散=扩散, 热源=热源),
        监测=''.join(_填(_监测项, 名=名, 面=面名, 操作=操作, 场=场, 权=权)
                     for 名, (面名, 操作, 场, 权) in 监测.items())
        + (''.join(_填(_截面监测项, 名=名, 点=_向量(点), 法向=_向量(垂直向量(流向)[0]))
                   for 名, 点 in 截面位置(流向, 多孔范围).items()) if 截面 else '')
        + (_yPlus if 模型 and 工况['壁面'] == '无滑移' and 'walls' in 面 else '')))
    写('system/fvSchemes', 'dictionary', _填(
        _fvSchemes, 湍流对流=''.join('    div(phi,%s)%s bounded Gauss upwind;\n' % (f, ' ' * (8 - len(f)))
                                    for f in 湍流场[模型]) if 模型 else '',
        壁面距离='wallDist        { method meshWave; }\n' if 模型 else ''))
    写('system/fvSolution', 'dictionary', _湍流求解设置(_fvSolution) if 模型 else _fvSolution)
    写('system/fvOptions', 'dictionary', 多孔)
    写('constant/transportProperties', 'dictionary',
       _填(_transportProperties, nu='%.9g' % 流体['有效运动粘度']))
    写('constant/turbulenceProperties', 'dictionary',
       _填(_turbulenceProperties_RAS, 模型=模型) if 模型 else _turbulenceProperties)
    T入 = '%.9g' % 工况['入口']['温度']
    if 模型:
        for 路径, 内容 in _湍流文件(湍流, 工况, 边界, T入).items():
            写(路径, 'volScalarField', 内容)
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
        'walls': 壁面温度条件(工况, 流体),
        'symmetry': 'type symmetry;'})))
    改 = ''.join(_填(_改类型, 面=名, 类=边界类型[名]) for 名 in ('symmetry', 'walls') if 名 in 面)
    编译 = 多孔范围 is not None and 工况['多孔']['换热']['方式'] == '水侧耦合'
    加层 = 要加边界层(工况, 面)
    if 加层:
        写('system/snappyHexMeshDict', 'dictionary', _填(
            _snappyHexMeshDict, 层数=边界层['层数'], 扩展比='%.9g' % 边界层['扩展比'],
            末层='%.9g' % 边界层['末层厚度比']))
    _写(os.path.join(算例目录, 'Allrun'), _填(_Allrun, 改类型=改, 加层=_加层命令 if 加层 else '',
                                            求解=_求解_编译 if 编译 else _求解_普通),
       可执行=True)
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
    if 名 in ('postProcessing', '结果.json', 'dynamicCode') or 名.startswith('log.'):
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


def generate_case(step路径, 边界, 网格, 输出目录, 工况=None, 步数=None, 流向=None):
    '''核对边界与网格一致后写出算例，返回 算例信息.json 的内容，另加 算例目录、算例信息路径 两项。

    不改动传入的字典、模型.step、标注.json 和第 4 步的 mesh.msh。
    输出目录里只替换本步的产物和旧的求解结果（log.*、postProcessing、时间步目录、结果.json）。
    '''
    步数 = _校验步数(步数)
    面标注, 多孔体 = _读边界(边界)
    if 工况 is None:
        if not isinstance(边界, dict) or '工况' not in 边界:
            raise 工况错误('没有工况: 边界不是第 2 步的返回值时，必须另外传入工况字典')
        工况 = 边界['工况']
    工况 = 校验工况(工况, bool(多孔体))
    if 流向 is None:
        流向 = 读流向(工况)
    n = 垂直向量(流向)[0]
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
    入口周长 = math.fsum(f.get('周长', 0.0) for f in 信息['面']
                         if f['编号'] in 外表面 and 面标注[f['编号']] == 'inlet')
    多孔范围 = None
    if 多孔集:                                   # 用多孔区各表面的中心估计它沿流向的起止位置
        s = [sum(a * b for a, b in zip(f['中心'], n)) for f in 信息['面'] if 多孔集 & set(f['相邻体'])]
        多孔范围 = (min(s), max(s))
    流体 = 流体参数(工况)
    湍流 = 湍流参数(工况, 流体, 入口面积, 入口周长)
    流量 = 工况['入口']['速度'] * 入口面积
    停留时间 = math.fsum(v['体积'] for v in 信息['体']) / 流量
    交界 = [f['面积'] for f in 信息['面']
            if f['编号'] in 内部面 and len(多孔集 & set(f['相邻体'])) == 1]
    流速, 长度, 压降, 出口温度 = None, None, None, None
    if len(交界) == 2:                           # 多孔区一进一出两个交界面时才能估算
        截面 = 0.5 * sum(交界)
        流速 = 流量 / 截面
        长度 = math.fsum(v['体积'] for v in 信息['体'] if v['编号'] in 多孔集) / 截面
        压降 = 多孔压降估算(工况, 流体, 流速, 长度)
        出口温度 = 出口温度估算(工况, 流体, 流速, 长度, 截面=截面)
    换 = 工况['多孔']['换热'] if 多孔集 else {}
    水侧 = 水侧摘要(换)
    if 水侧 is not None:
        if 多孔范围 is None or not 多孔范围[1] > 多孔范围[0]:
            raise 算例错误('水侧耦合要求多孔区沿流向有长度')
        水侧['水出口温度估算'] = 水侧估算(工况, 流体, 流速, 长度, 截面)[1] if 压降 is not None else None
    摘要 = {'边界': {名: 边界类型[名] for 名 in sorted(面名)}, '多孔区': 多孔区名 if 多孔集 else None,
            '步数': 步数, '时间步': 伪时间步(工况, 流体, 停留时间), '停留时间': 停留时间, '流向': n, '监测量': sorted(list(监测) + (list(截面监测) if 压降 is not None else [])),
            '入口面积': 入口面积, '入口速度': 工况['入口']['速度'], '入口温度': 工况['入口']['温度'],
            '出口表压': 工况['出口']['表压'], '密度': 流体['密度'], '比热': 流体['比热'],
            '分子粘度': 流体['分子粘度'], '有效运动粘度': 流体['有效运动粘度'],
            '密度模型': 流体['密度模型'], '粘度模型': 流体['粘度模型'],
            '有效导热系数': 流体['有效导热系数'], '湍流': 湍流,
            '边界层': dict(边界层) if 要加边界层(工况, 面名) else None,
            '壁面换热': 工况.get('壁面换热', {'方式': '绝热'}),
            '理论入口体积流量': 流量, '理论入口质量流量': 流量 * 流体['密度'],
            '多孔范围': list(多孔范围) if 多孔范围 else None, '多孔区流速': 流速, '多孔区长度': 长度,
            '多孔压降估算': 压降, '出口温度估算': 出口温度, '水侧': 水侧,
            '压力截面': ({名: list(点) for 名, 点 in 截面位置(n, 多孔范围).items()}
                         if 压降 is not None else None)}

    输出目录 = os.fspath(输出目录) or '.'
    os.makedirs(输出目录, exist_ok=True)
    临时 = os.path.join(输出目录, '.case.' + uuid.uuid4().hex + '.tmp')
    try:
        写算例(临时, 工况, 面名, 多孔范围, 流向, 步数, 截面=压降 is not None, 停留时间=停留时间, 湍流=湍流)
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
