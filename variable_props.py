# -*- coding: utf-8 -*-
'''变物性（密度、粘度、导热系数随温度变化）+ 重力（浮升力）的算例，见 docs/12 的“变物性与浮升力”。

工况里写 '物性模型': '变物性' 时由 generate_case 调用，求解器换成 buoyantSimpleFoam（可压缩定常 SIMPLE，
解焓 h，p = p_rgh + ρ·g·x）。不写 物性模型 或写 '常物性' 时完全不经过这里，原来的算例不变。

物性:
    密度    理想气体（工作压力 + 摩尔质量）: ρ = p工作·M/(R·T)，OpenFOAM 用 incompressiblePerfectGas
            （密度只随温度变，不随当地压力变，低马赫数流动的通常做法）；
            常数 密度: 不随温度变（icoPolynomial 的常数项）；温度表 的 密度 列: 按表拟合多项式（icoPolynomial）。
    粘度    Sutherland（默认空气系数）: 在算例温度范围上拟合成 4 次多项式（OpenFOAM 的 polynomial 输运）；
            常数 分子粘度；温度表 的 分子粘度 列: 拟合多项式。
    导热    常数 导热系数；分子普朗特数（气体默认 0.7）: κ = μ(T)·cp/Pr，与 μ 同一个多项式乘常数；
            温度表 的 导热系数 列: 拟合多项式。
    比热    常数（hPolynomial 的常数项），h = cp·(T − 298.15)。
温度表的多项式阶数取 min(点数 − 1, 4)；点数不超过 5 时就是过各点的插值多项式。
拟合误差（相对）记进 算例信息.json 的 物性拟合。
重力: 工况的 重力 [gx, gy, gz]，默认 [0, 0, -9.81]（z 轴向上）；[0, 0, 0] 时没有浮升力，只有变物性。
'''
import math

from generate_boundary import (工况错误, 默认气体普朗特数, 默认湍流普朗特数, 体积换热系数, 壁温两端, 湍流模型名)

R气体 = 8.314462
标准大气压 = 101325.0
Sutherland阶数 = 4
_采样点数 = 41
表阶数 = 4
拟合误差限 = 0.02          # 温度表的多项式在表里各点的最大相对误差超过它时报工况错误
求解器 = 'buoyantSimpleFoam'
日志名 = 'log.' + 求解器


# ---------------------------------------------------------------- 多项式拟合
def _解线性方程(A, b):
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        主 = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[主] = M[主], M[c]
        if M[c][c] == 0:
            raise 工况错误('物性多项式拟合失败（方程奇异）')
        for r in range(n):
            if r != c:
                k = M[r][c] / M[c][c]
                for j in range(c, n + 1):
                    M[r][j] -= k * M[c][j]
    return [M[i][n] / M[i][i] for i in range(n)]


def 拟合(温度, 值, 阶):
    '''按相对误差最小二乘拟合 值 ≈ Σ a_i·T^i（i = 0..阶），返回 [a_0, ..., a_阶]。
    先在缩放变量 x = (T − Tc)/Th ∈ [−1, 1] 上拟合（避免病态），再展开成 T 的幂。'''
    阶 = min(阶, len(温度) - 1)
    Tc = 0.5 * (min(温度) + max(温度))
    Th = 0.5 * (max(温度) - min(温度)) or 1.0
    x = [(t - Tc) / Th for t in 温度]
    w = [1.0 / v for v in 值]                  # 相对误差
    A = [[math.fsum(w[k] ** 2 * x[k] ** (i + j) for k in range(len(x))) for j in range(阶 + 1)]
         for i in range(阶 + 1)]
    b = [math.fsum(w[k] ** 2 * x[k] ** i * 值[k] for k in range(len(x))) for i in range(阶 + 1)]
    c = _解线性方程(A, b)
    a = [0.0] * (阶 + 1)
    for j, cj in enumerate(c):                 # cj·((T − Tc)/Th)^j = cj/Th^j · Σ C(j,i)·T^i·(−Tc)^(j−i)
        for i in range(j + 1):
            a[i] += cj / Th ** j * math.comb(j, i) * (-Tc) ** (j - i)
    return a


def 多项式值(a, T):
    值 = 0.0
    for ai in reversed(a):
        值 = 值 * T + ai
    return 值


def _八项(a):
    '''OpenFOAM Polynomial<8> 的 8 个系数。'''
    a = list(a) + [0.0] * (8 - len(a))
    return '(' + ' '.join('%.17g' % (x + 0.0) for x in a) + ')'


def sutherland(T, 参考粘度=1.716e-5, 参考温度=273.15, 常数=110.4):
    return 参考粘度 * (T / 参考温度) ** 1.5 * (参考温度 + 常数) / (T + 常数)


# ---------------------------------------------------------------- 物性模型
def 工况温度(工况):
    '''工况里给定的温度（入口、壁温、多孔区壁温、水入口温度），K。'''
    温度 = [工况['入口']['温度']]
    换 = 工况.get('壁面换热', {})
    if 换.get('方式') == '定壁温':
        温度.append(换['壁温'])
    孔换 = (工况.get('多孔') or {}).get('换热', {})
    if 孔换.get('方式') == '温差':
        温度 += list(壁温两端(孔换))
    elif 孔换.get('方式') == '水侧耦合':
        温度.append(孔换['水入口温度'])
    return 温度


def 物性(工况, 温度):
    '''变物性的多项式和 Python 里用的物性函数。温度: 算例里会出现的温度（工况温度 + 一维估算的出口温度）。
    返回 {'密度': f(T), '粘度': f(T), '导热': f(T), '比热': cp, 'OpenFOAM': {...}, '拟合': {...}}。'''
    物 = 工况['物性']
    表 = 物.get('温度表')
    T低, T高 = min(温度), max(温度)
    if 表:
        表T = 表['温度']
        if T低 < 表T[0] - 1e-9 or T高 > 表T[-1] + 1e-9:
            raise 工况错误('物性.温度表 的温度范围 %g–%g K 没有覆盖算例的温度 %g–%g K（入口、壁温、水温和一维估算的'
                           '出口温度），请把表补到这个范围' % (表T[0], 表T[-1], T低, T高))
    范围 = (0.8 * T低, 1.3 * T高)          # Sutherland 拟合的温度范围，两头留余量给局部的高低温
    拟合信息 = {}
    cp = 物['比热']

    def 按表(列):
        a = 拟合(表['温度'], 表[列], 表阶数)
        误差 = max(abs(多项式值(a, t) / v - 1) for t, v in zip(表['温度'], 表[列]))
        if 误差 > 拟合误差限:
            raise 工况错误('物性.温度表 的 %s 列用 %d 次多项式拟合，表里的点最大相对误差 %.2g%%，超过 %g%%；'
                           '表里的值可能有错或变化太剧烈' % (列, len(a) - 1, 误差 * 100, 拟合误差限 * 100))
        样 = [表['温度'][0] + (表['温度'][-1] - 表['温度'][0]) * i / 200 for i in range(201)]
        if min(多项式值(a, t) for t in 样) <= 0:
            raise 工况错误('物性.温度表 的 %s 列拟合出的多项式在表的温度范围内出现非正值，请检查表' % 列)
        拟合信息[列] = {'来源': '温度表', '阶数': len(a) - 1, '温度范围': [表['温度'][0], 表['温度'][-1]],
                        '最大相对误差': 误差}
        return a

    of = {}
    if 表 and '密度' in 表:
        a = 按表('密度')
        of['状态方程'] = ('icoPolynomial', 'rhoCoeffs<8> %s;' % _八项(a))
        密度 = lambda T, a=a: 多项式值(a, T)
        拟合信息['密度']['模型'] = '温度表多项式'
    elif '密度' in 物:
        ρ = 物['密度']
        of['状态方程'] = ('icoPolynomial', 'rhoCoeffs<8> %s;' % _八项([ρ]))
        密度 = lambda T, ρ=ρ: ρ
        拟合信息['密度'] = {'模型': '常数', '值': ρ}
    else:
        p, M = 物['工作压力'], 物['摩尔质量']
        of['状态方程'] = ('incompressiblePerfectGas', 'pRef %.12g;' % p)
        密度 = lambda T, p=p, M=M: p * M / (R气体 * T)
        拟合信息['密度'] = {'模型': '理想气体 ρ = p工作·M/(R·T)', '工作压力': p, '摩尔质量': M}
    of['摩尔质量'] = 物.get('摩尔质量', 0.018015) * 1000      # g/mol；液体只用于 OpenFOAM 的 R，不影响结果

    if 表 and '分子粘度' in 表:
        a粘 = 按表('分子粘度')
        拟合信息['分子粘度']['模型'] = '温度表多项式'
    elif '分子粘度' in 物:
        a粘 = [物['分子粘度']]
        拟合信息['分子粘度'] = {'模型': '常数', '值': 物['分子粘度']}
    else:
        系数 = (物.get('萨瑟兰参考粘度', 1.716e-5), 物.get('萨瑟兰参考温度', 273.15), 物.get('萨瑟兰常数', 110.4))
        样 = [0.5 * (范围[0] + 范围[1]) - 0.5 * (范围[1] - 范围[0]) * math.cos(math.pi * (k + 0.5) / _采样点数)
             for k in range(_采样点数)]
        a粘 = 拟合(样, [sutherland(t, *系数) for t in 样], Sutherland阶数)
        检 = [范围[0] + (范围[1] - 范围[0]) * i / 400 for i in range(401)]
        拟合信息['分子粘度'] = {'模型': 'Sutherland（%d 次多项式拟合）' % Sutherland阶数, '温度范围': list(范围),
                              '最大相对误差': max(abs(多项式值(a粘, t) / sutherland(t, *系数) - 1) for t in 检)}

    if 表 and '导热系数' in 表:
        a导 = 按表('导热系数')
        拟合信息['导热系数']['模型'] = '温度表多项式'
    elif '导热系数' in 物:
        a导 = [物['导热系数']]
        拟合信息['导热系数'] = {'模型': '常数', '值': 物['导热系数']}
    else:
        Pr = 物.get('分子普朗特数', 默认气体普朗特数)
        a导 = [x * cp / Pr for x in a粘]
        拟合信息['导热系数'] = {'模型': 'μ(T)·cp/Pr', '分子普朗特数': Pr}
    of['muCoeffs'] = _八项(a粘)
    of['kappaCoeffs'] = _八项(a导)
    of['CpCoeffs'] = _八项([cp])
    return {'密度': 密度, '粘度': lambda T: 多项式值(a粘, T), '导热': lambda T: 多项式值(a导, T), '比热': cp,
            'OpenFOAM': of, '拟合': 拟合信息, '温度范围': [T低, T高]}


def 流体参数(工况, 模):
    '''入口温度下的物性，键与 generate_case.流体参数 相同（变物性没有常数湍流粘度）。'''
    T = 工况['入口']['温度']
    ρ, μ, κ, cp = 模['密度'](T), 模['粘度'](T), 模['导热'](T), 模['比热']
    return {'密度': ρ, '分子粘度': μ, '有效运动粘度': μ / ρ, '比热': cp, '热扩散系数': κ / (ρ * cp),
            '有效导热系数': κ, '分子导热系数': κ, '密度模型': 模['拟合']['密度']['模型'],
            '粘度模型': 模['拟合']['分子粘度']['模型']}


def 参考压力(工况):
    '''p 场的绝对压力零点（Pa）: 理想气体是 工作压力，其余取 101325。出口绝对压力 = 参考压力 + 出口表压。'''
    return 工况['物性'].get('工作压力', 标准大气压)


# ---------------------------------------------------------------- 一维估算
def 温度剖面(工况, G, 长度, 截面, cp, 细分=2000):
    '''均匀流过多孔区时沿流向 细分+1 个节点的烟气温度。G: 质量流量密度 kg/(m²·s)（不随位置变）。
    能量方程只含 G·cp，与密度无关，所以出口温度和常物性的一维估算相同。水侧耦合按 ε-NTU 的水出口温度
    定出 s = 0 处的水温，再把烟气和水的两个方程一起往下游推（RK4）。'''
    换, T = 工况['多孔']['换热'], 工况['入口']['温度']
    ds = 长度 / 细分
    G热容 = G * cp
    剖面 = [T]
    if 换['方式'] == '水侧耦合':
        from generate_case import 水侧估算
        流体 = {'密度': 1.0, '比热': cp}
        Tw出 = 水侧估算(工况, 流体, G, 长度, 截面)[1]
        逆流 = 换['水流动'] == '逆流'
        Tw = Tw出 if 逆流 else 换['水入口温度']
        kg = 体积换热系数(换) / G热容
        kw = 体积换热系数(换) * 截面 / (换['水质量流量'] * 换['水比热']) * (-1 if 逆流 else 1)

        def f(y):
            q = y[0] - y[1]
            return (-kg * q, kw * q)
        y = (T, Tw)
        for _ in range(细分):
            k1 = f(y)
            k2 = f((y[0] + 0.5 * ds * k1[0], y[1] + 0.5 * ds * k1[1]))
            k3 = f((y[0] + 0.5 * ds * k2[0], y[1] + 0.5 * ds * k2[1]))
            k4 = f((y[0] + ds * k3[0], y[1] + ds * k3[1]))
            y = (y[0] + ds / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
                 y[1] + ds / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))
            剖面.append(y[0])
        return 剖面
    for n in range(细分):
        if 换['方式'] == '温差':
            T入侧, T出侧 = 壁温两端(换)
            Tw = T入侧 + (T出侧 - T入侧) * (n + 0.5) / 细分
            a = 换['h_av'] * ds / G热容
            T = (T * (1 - a / 2) + a * Tw) / (1 + a / 2)
        else:
            T = T + 换['体积热源'] * ds / G热容
        剖面.append(T)
    return 剖面


def 多孔压降(工况, 模, G, 长度, 截面, 细分=2000):
    '''均匀流过多孔区的流动压降（Pa，不含重力），返回 (压降, 加速压降, 出口温度, 进口密度, 出口密度)。
    dp/ds = −(d·μ(T)·u + f·½ρ(T)·u²) − d(G·u)/ds，u = G/ρ(T)；
    加速项积分为 G²·(1/ρ出 − 1/ρ进)，是截面压差里由密度变化引起的那部分（气体冷却时为负）。'''
    孔 = 工况['多孔']
    剖面 = 温度剖面(工况, G, 长度, 截面, 模['比热'], 细分)
    ds = 长度 / 细分
    阻力 = 0.0
    for a, b in zip(剖面, 剖面[1:]):
        T = 0.5 * (a + b)
        ρ = 模['密度'](T)
        u = G / ρ
        阻力 += ds * (孔['粘性阻力系数'] * 模['粘度'](T) * u + 孔['惯性阻力系数'] * 0.5 * ρ * u * u)
    ρ进, ρ出 = 模['密度'](剖面[0]), 模['密度'](剖面[-1])
    加速 = G * G * (1 / ρ出 - 1 / ρ进)
    return 阻力 + 加速, 加速, 剖面[-1], ρ进, ρ出


# ---------------------------------------------------------------- 算例文件
_thermophysicalProperties = '''thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       polynomial;
    thermo          hPolynomial;
    equationOfState @状态方程@;
    specie          specie;
    energy          sensibleEnthalpy;
}

mixture
{
    specie
    {
        molWeight       @摩尔质量@;
    }
    equationOfState
    {
        @状态参数@
    }
    thermodynamics
    {
        Hf              0;
        Sf              0;
        CpCoeffs<8>     @Cp@;
    }
    transport
    {
        muCoeffs<8>     @mu@;
        kappaCoeffs<8>  @kappa@;
    }
}
'''
_g = '''dimensions      [0 1 -2 0 0 0 0];
value           @g@;
'''
# p 是绝对压力（约 1e5 Pa），监测文件要 12 位有效数字才能分辨 0.001 Pa 以下的压差
_controlDict = '''application     buoyantSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         @步数@;
deltaT          1;
writeControl    timeStep;
writeInterval   @步数@;
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;

functions
{
@监测@}
'''
# 梯度用 leastSquares：buoyantSimpleFoam 的动量方程用面法向压力梯度重构（fvc::reconstruct），
# 在四面体网格上配 Gauss linear 时，多孔区压降稳定偏低约 5.6%（与网格疏密、松弛、
# 非正交修正都无关）；改用 leastSquares 后与一维估算差 0.1%。见 docs/12 第 10 节。
_fvSchemes = '''ddtSchemes      { default steadyState; }
gradSchemes
{
    default         leastSquares;
    grad(U)         cellLimited Gauss linear 1;
}
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;
@湍流对流@    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
wallDist        { method meshWave; }
'''
# 松弛因子取 OpenFOAM 教程 buoyantCavity 的值（p_rgh 0.7，U、h 0.3）。按不可压的 p 0.3 / U 0.7 时，
# 700 K 烟气在多孔区里冷却的算例第 3 步就出现负温度发散（docs/12）。
_fvSolution = '''solvers
{
    p_rgh
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    "(U|k|epsilon|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.1;
    }
    h
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-10;
        relTol          0.001;
    }
}
SIMPLE
{
    momentumPredictor yes;
    nNonOrthogonalCorrectors 2;
    residualControl
    {
        p_rgh           1e-6;
        U               1e-6;
        h               1e-8;
        "(k|epsilon|omega)" 1e-6;
    }
}
relaxationFactors
{
    fields      { rho 1; p_rgh 0.7; }
    equations   { "(U|k|epsilon|omega)" 0.3; h 0.3; }
}
'''
_turbulenceProperties_RAS = '''simulationType  RAS;

RAS
{
    RASModel        @模型@;
    turbulence      on;
    printCoeffs     on;
    @模型@Coeffs
    {
        Prt         @Prt@;
    }
}
'''
_定热源 = '''heat
{
    type            scalarSemiImplicitSource;
    selectionMode   cellZone;
    cellZone        porous_zone;
    volumeMode      specific;
    sources
    {
        h           (@Su@ 0);
    }
}
'''
# 温差换热: 源项 h_av·(Tw − T)（W/m³）写成焓的形式 kh·(hw − h)，kh = h_av/cp，hw = h + cp·(Tw − T)，
# −kh·h 放进对角线（真正隐式）。T 是上一次迭代的温度，与 h 一一对应，收敛后与温差形式完全相同。
_温差源 = '''wallTemperature
{
    type            scalarCodedSource;
    name            wallTemperature;
    selectionMode   cellZone;
    cellZone        porous_zone;
    fields          (h);

    codeInclude
    #{
    #};

    codeCorrect
    #{
    #};

    codeConstrain
    #{
    #};

    codeAddSup
    #{
    #};

    codeAddSupRho
    #{
        const scalar kh = @kh@;
        const scalar Cp = @Cp@;
        const scalar T0 = @T入侧@;
        const scalar T1 = @T出侧@;
        const vector n@n@;
        const scalar s0 = @起@;
        const scalar L = @长@;
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        const scalarField& h = eqn.psi().primitiveField();
        const scalarField& V = mesh_.V();
        const vectorField& C = mesh_.C();
        scalarField& src = eqn.source();
        scalarField& dg = eqn.diag();
        forAll(cells_, i)
        {
            const label c = cells_[i];
            const scalar s = min(max(((C[c] & n) - s0)/L, scalar(0)), scalar(1));
            const scalar Tw = T0 + (T1 - T0)*s;
            src[c] -= kh*(h[c] + Cp*(Tw - T[c]))*V[c];
            dg[c] -= kh*V[c];
        }
    #};
}
'''


def _水侧源(工况, 流体, 流向, 多孔范围):
    '''常物性的水侧 scalarCodedSource 改成焓方程: 作用在 h 上，代码放进 codeAddSupRho（buoyantSimpleFoam 调的是
    带 ρ 的 addSup），温度从 T 场取，源项 kh·(hw − h)，kh = h_av/cp，hw = h + cp·(Tw − T)。'''
    import generate_case as gc
    文本 = gc.水侧源项(工况, 流体, 流向, 多孔范围)
    cp = 流体['比热']
    hav = 体积换热系数(工况['多孔']['换热'])
    for 旧, 新 in (
            ('fields          (T);', 'fields          (h);'),
            ('                codeAddSup\n                #{\n',
             '                codeAddSup\n                #{\n                #};\n\n'
             '                codeAddSupRho\n                #{\n'),
            ('const scalar k = %.12g;' % (hav / (流体['密度'] * cp)),
             'const scalar kh = %.12g;\n                    const scalar Cp = %.12g;' % (hav / cp, cp)),
            ('const scalarField& T = eqn.psi().primitiveField();',
             'const volScalarField& T = mesh_.lookupObject<volScalarField>("T");\n'
             '                    const scalarField& hh = eqn.psi().primitiveField();'),
            ('src[cs[i]] -= k*Tw*V[cs[i]];\n                        dg[cs[i]] -= k*V[cs[i]];',
             'src[cs[i]] -= kh*(hh[cs[i]] + Cp*(Tw - T[cs[i]]))*V[cs[i]];\n'
             '                        dg[cs[i]] -= kh*V[cs[i]];')):
        assert 文本.count(旧) == 1, 旧
        文本 = 文本.replace(旧, 新)
    return 文本


def _缩进(文本, 空格):
    return ''.join((' ' * 空格 + 行 if 行.strip() else 行) for 行 in 文本.splitlines(True))


def 壁面温度条件(工况):
    换 = 工况.get('壁面换热', {'方式': '绝热'})
    if 换['方式'] == '定壁温':
        return 'type fixedValue; value uniform %.9g;' % 换['壁温']
    if 换['方式'] == '定热流':     # 按热物性（含湍流 alphat）把热流密度换成壁面温度梯度
        return ('type externalWallHeatFluxTemperature; mode flux; q uniform %.9g; kappaMethod fluidThermo; '
                'value uniform %.9g;' % (换['热流密度'], 工况['入口']['温度']))
    return 'type zeroGradient;'


_限温 = '''limitT
{
    type            limitTemperature;
    active          yes;
    selectionMode   all;
    min             @低@;
    max             @高@;
}
'''


def 温度限(模):
    '''求解中温度的上下限（limitTemperature）: 迭代初期温度可能远超物理范围，多项式外推出负密度会让求解崩溃。
    按温度表时就是表的范围（表外的多项式外推不可信）；否则取算例温度范围的 0.5–2 倍。'''
    表范围 = [v['温度范围'] for v in 模['拟合'].values() if v.get('来源') == '温度表']
    if 表范围:
        return max(r[0] for r in 表范围), min(r[1] for r in 表范围)
    T低, T高 = 模['温度范围']
    return 0.5 * T低, 2.0 * T高


def 写算例(算例目录, 工况, 面, 多孔范围, 流向, 步数, 截面, 湍流, 模, 流体):
    '''写 buoyantSimpleFoam 算例的 0/、constant/、system/ 和 Allrun。参数同 generate_case.写算例，
    另加 模（物性() 的返回值）和 流体（入口物性）。'''
    import generate_case as gc
    面 = set(面)
    模型 = 湍流模型名(工况)
    p0 = 参考压力(工况)
    p出 = p0 + 工况['出口']['表压']
    无滑移 = 工况['壁面'] == '无滑移'

    def 写(相对路径, 类, 内容):
        import os
        路径 = os.path.join(算例目录, *相对路径.split('/'))
        gc._写(路径, gc._填(gc._头, 类=类, 名=os.path.basename(路径)) + 内容)

    def 边界(条件):
        return ''.join(gc._填(gc._块, 面=名, 内容=条件[名]) for 名 in sorted(面))

    fvOptions = ''
    编译 = False
    if 多孔范围 is not None:
        e1, e2 = gc.垂直向量(流向)
        d, f = gc.阻力系数(工况, 流体)           # 变物性没有常数湍流粘度，d 不折算
        fvOptions = gc._填(gc._多孔, d=gc._向量(d), f=gc._向量(f), e1=gc._向量(e1), e2=gc._向量(e2))
        换 = 工况['多孔']['换热']
        起, 终 = 多孔范围
        if not 终 > 起:
            raise gc.算例错误('多孔区沿流向的长度必须大于 0')
        if 换['方式'] == '定热源':
            fvOptions += '\n' + gc._填(_定热源, Su='%.12g' % 换['体积热源'])
        elif 换['方式'] == '温差':
            T入侧, T出侧 = 壁温两端(换)
            fvOptions += '\n' + gc._填(
                _温差源, kh='%.12g' % (换['h_av'] / 流体['比热']), Cp='%.12g' % 流体['比热'],
                T入侧='%.12g' % T入侧, T出侧='%.12g' % T出侧,
                n='(%s)' % ', '.join('%.12g' % (x + 0.0) for x in gc.垂直向量(流向)[0]),
                起='%.12g' % 起, 长='%.12g' % (终 - 起))
            编译 = True
        else:
            水 = _水侧源(工况, 流体, 流向, 多孔范围)
            fvOptions += '\n' + _缩进(''.join(行[12:] if 行.startswith(' ' * 12) else 行
                                            for 行 in 水.splitlines(True)), 0)
            编译 = True

    T限 = 温度限(模)
    fvOptions += ('\n' if fvOptions else '') + gc._填(_限温, 低='%.9g' % T限[0], 高='%.9g' % T限[1])
    n = gc.垂直向量(流向)[0]
    监测项 = dict(监测)
    写('system/controlDict', 'dictionary', gc._填(
        _controlDict, 步数=步数,
        监测=''.join(gc._填(gc._监测项, 名=名, 面=面名, 操作=操作, 场=场, 权=权)
                     for 名, (面名, 操作, 场, 权) in 监测项.items())
        + (''.join(gc._填(gc._截面监测项, 名=名, 点=gc._向量(点), 法向=gc._向量(n)).replace(
            'fields          (p);', 'fields          (%s);' % 场)
            for 名, 场, 点 in 截面监测项(n, 多孔范围)) if 截面 else '')
        + (gc._yPlus if 模型 and 无滑移 and 'walls' in 面 else '')))
    写('system/fvSchemes', 'dictionary', gc._填(
        _fvSchemes, 湍流对流=''.join('    div(phi,%s)%s bounded Gauss upwind;\n' % (f, ' ' * (8 - len(f)))
                                    for f in gc.湍流场[模型]) if 模型 else ''))
    写('system/fvSolution', 'dictionary', _fvSolution)
    写('system/fvOptions', 'dictionary', fvOptions)
    of = 模['OpenFOAM']
    写('constant/thermophysicalProperties', 'dictionary', gc._填(
        _thermophysicalProperties, 状态方程=of['状态方程'][0], 状态参数=of['状态方程'][1],
        摩尔质量='%.12g' % of['摩尔质量'], Cp=of['CpCoeffs'], mu=of['muCoeffs'], kappa=of['kappaCoeffs']))
    写('constant/g', 'uniformDimensionedVectorField', gc._填(_g, g=gc._向量(工况['重力'])))
    写('constant/turbulenceProperties', 'dictionary',
       gc._填(_turbulenceProperties_RAS, 模型=模型, Prt='%.9g' % 工况['物性'].get('湍流普朗特数', 默认湍流普朗特数))
       if 模型 else gc._turbulenceProperties)
    T入 = '%.9g' % 工况['入口']['温度']
    if 模型:
        for 路径, 内容 in gc._湍流文件(湍流, 工况, 边界, T入).items():
            写(路径, 'volScalarField', 内容)
        Prt = '%.9g' % 工况['物性'].get('湍流普朗特数', 默认湍流普朗特数)
        写('0/alphat', 'volScalarField', gc._填(gc._场, 量纲='[1 -1 -1 0 0 0 0]', 初值='0', 边界=边界({
            'inlet': 'type calculated; value uniform 0;',
            'outlet': 'type calculated; value uniform 0;',
            'walls': ('type compressible::alphatWallFunction; Prt %s; value uniform 0;' % Prt
                      if 无滑移 else 'type calculated; value uniform 0;'),
            'symmetry': 'type symmetry;'})))
    U初 = gc._向量([工况['入口']['速度'] * x for x in 流向])       # 初场取入口速度: 变物性时 U=0 起步温度会失控
    写('0/U', 'volVectorField', gc._填(gc._场, 量纲='[0 1 -1 0 0 0 0]', 初值=U初, 边界=边界({
        'inlet': 'type surfaceNormalFixedValue; refValue uniform %.9g;' % -工况['入口']['速度'],
        'outlet': 'type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0);',
        'walls': 'type slip;' if 工况['壁面'] == '滑移' else 'type noSlip;',
        'symmetry': 'type symmetry;'})))
    P = '%.12g' % p出
    写('0/p', 'volScalarField', gc._填(gc._场, 量纲='[1 -1 -2 0 0 0 0]', 初值=P, 边界=边界({
        'inlet': 'type calculated; value uniform %s;' % P,
        'outlet': 'type calculated; value uniform %s;' % P,
        'walls': 'type calculated; value uniform %s;' % P,
        'symmetry': 'type symmetry;'})))
    写('0/p_rgh', 'volScalarField', gc._填(gc._场, 量纲='[1 -1 -2 0 0 0 0]', 初值=P, 边界=边界({
        'inlet': 'type fixedFluxPressure; value uniform %s;' % P,
        'outlet': 'type prghPressure; p uniform %s; value uniform %s;' % (P, P),
        'walls': 'type fixedFluxPressure; value uniform %s;' % P,
        'symmetry': 'type symmetry;'})))
    写('0/T', 'volScalarField', gc._填(gc._场, 量纲='[0 0 0 1 0 0 0]', 初值=T入, 边界=边界({
        'inlet': 'type fixedValue; value uniform %s;' % T入,
        'outlet': 'type inletOutlet; inletValue uniform %s; value uniform %s;' % (T入, T入),
        'walls': 壁面温度条件(工况),
        'symmetry': 'type symmetry;'})))
    改 = ''.join(gc._填(gc._改类型, 面=名, 类=gc.边界类型[名]) for 名 in ('symmetry', 'walls') if 名 in 面)
    加层 = gc.要加边界层(工况, 面)
    if 加层:
        写('system/snappyHexMeshDict', 'dictionary', gc._填(
            gc._snappyHexMeshDict, 层数=gc.边界层['层数'], 扩展比='%.9g' % gc.边界层['扩展比'],
            末层='%.9g' % gc.边界层['末层厚度比']))
    求解 = (gc._求解_编译 if 编译 else gc._求解_普通).replace('simpleFoam > log.simpleFoam',
                                                         '%s > %s' % (求解器, 日志名))
    import os
    gc._写(os.path.join(算例目录, 'Allrun'), gc._填(gc._Allrun, 改类型=改, 加层=gc._加层命令 if 加层 else '',
                                                  求解=求解), 可执行=True)
    return 流体


# 变物性的监测量: p 是绝对压力（Pa），phi 是质量流量（kg/s），另加进出口和截面的平均密度
监测 = {'inletP': ('inlet', 'areaAverage', 'p', ''), 'inletFlow': ('inlet', 'sum', 'phi', ''),
        'outletFlow': ('outlet', 'sum', 'phi', ''),
        'outletT': ('outlet', 'weightedAverage', 'T', 'weightField     phi;'),
        'inletRho': ('inlet', 'areaAverage', 'rho', ''), 'outletRho': ('outlet', 'areaAverage', 'rho', '')}
截面监测 = {'porousInP': ('起', 'p'), 'porousOutP': ('终', 'p'), 'porousInRho': ('起', 'rho'),
            'porousOutRho': ('终', 'rho')}


def 截面监测项(n, 多孔范围):
    import generate_case as gc
    点 = {gc.截面监测[名]: 位置 for 名, 位置 in gc.截面位置(n, 多孔范围).items()}
    return [(名, 场, 点[端]) for 名, (端, 场) in 截面监测.items()]
