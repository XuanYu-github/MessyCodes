import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math

# 设置全局字体为支持中文的字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 使用 SimHei 字体（黑体）
matplotlib.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# --- 参数定义 ---
# 基本常数
beta = 0.006019
Lambda = 0.00002  # s
lambda_ = 0.150  # s^-1
lambda_I = 2.9e-5  # s^-1
lambda_X = 2.1e-5  # s^-1
gamma_I = 0.059
gamma_X = 0.023
sigma_aX = 3.5e-18  # cm^2

# 反应堆参数 (来自表 1)
P0 = 3000.0  # MW
Eff = 3.2e-11  # MWs 
D = 316.0  # cm
h = 355.0  # cm

# --- 计算衍生参数 ---
# 计算堆芯体积 V
V = 2 * math.pi * (D / 2)**2 * h  # cm^3

# 计算初始平均中子注量率 phi_0
phi_0 = P0 / (Eff * V)  # cm^-2 s^-1

# 计算 theta 和 vartheta (根据公式 3)
theta = gamma_X / gamma_I
vartheta = (sigma_aX * phi_0) / lambda_X

# 打印计算出的关键参数以供检查
print(f"计算得到的参数值:")
print(f"  V = {V:.4e} cm^3")
print(f"  phi_0 = {phi_0:.4e} cm^-2 s^-1")
print(f"  theta = {theta:.4f}")
print(f"  vartheta = {vartheta:.4e}") # 使用科学计数法显示更精确的小值

# --- 反应性函数 ---
def rho_t(t, case):
    """根据时间和情景计算反应性 rho(t)"""
    if case == 1:  # 正向阶跃
        return 0 if t < 10 else 0.001
    elif case == 2:  # 负向阶跃
        return 0 if t < 10 else -0.001
    elif case == 3:  # 正弦变化
        return 0 if t < 10 else -0.001 + 0.0002 * np.sin(t - 10)
    else:
        return 0 # 默认或错误情况

# --- 微分方程组定义 ---
def derivatives(t, y, case):
    """计算状态向量 y 在时间 t 的导数 dy/dt
    y = [n_r, c_r, I_r, X_r]
    """
    n_r, c_r, I_r, X_r = y
    # 添加检查，防止 inf/nan 传播
    if not np.all(np.isfinite(y)):
        # 如果输入已经是无效值，返回零向量或 NaN 向量，避免进一步计算和警告
        # 返回零向量可能使曲线在此后变平，返回 NaN 会使曲线中断
        # return np.zeros_like(y)
        return np.full_like(y, np.nan)

    rho = rho_t(t, case)

    # 点堆动力学方程 (1)
    dn_r_dt = ((rho - beta) / Lambda) * n_r + (beta / Lambda) * c_r
    dc_r_dt = lambda_ * (n_r - c_r)

    # 氙碘方程 (2)
    dI_r_dt = lambda_I * (n_r - I_r)
    term1 = (I_r + theta * n_r) / (1 + theta)
    term2 = (1 + vartheta * n_r) / (1 + theta) * X_r
    dX_r_dt = lambda_X * (1 + vartheta) * (term1 - term2)

    result = np.array([dn_r_dt, dc_r_dt, dI_r_dt, dX_r_dt])
    # 再次检查计算结果是否有效
    if not np.all(np.isfinite(result)):
         print(f"Warning: Invalid derivative calculated at t={t}, y={y}") # 调试信息
         return np.full_like(y, np.nan) # 如果计算结果无效，返回 NaN
    return result

# --- 四阶龙格-库塔步进函数 ---
def rk4_step(t, y, dt, case):
    """执行一个 RK4 步进"""
    k1 = derivatives(t, y, case)
    # 检查 k1 是否有效
    if np.any(np.isnan(k1)): return y + np.nan

    k2 = derivatives(t + dt/2, y + dt/2 * k1, case)
    if np.any(np.isnan(k2)): return y + np.nan

    k3 = derivatives(t + dt/2, y + dt/2 * k2, case)
    if np.any(np.isnan(k3)): return y + np.nan

    k4 = derivatives(t + dt, y + dt * k3, case)
    if np.any(np.isnan(k4)): return y + np.nan

    y_next = y + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)
    return y_next

# --- 求解函数 ---
def solve(case, n_r0, t_max=400, dt=0.1): 
    """使用 RK4 求解给定初始条件和情景下的 ODE"""
    print(f"  Solving case {case} for n_r0={n_r0} with dt={dt}...") # 添加打印信息
    t = np.arange(0, t_max + dt, dt)
    y = np.zeros((len(t), 4)) # 存储 [n_r, c_r, I_r, X_r]

    # 设置初始条件 (基于稳态假设 rho(0)=0)
    c_r0 = n_r0
    I_r0 = n_r0
    if (1 + vartheta * n_r0) == 0:
         X_r0 = 0
    else:
        X_r0 = n_r0 * (1 + theta) / (1 + vartheta * n_r0)
    y[0] = [n_r0, c_r0, I_r0, X_r0]

    # 时间步进求解
    for i in range(len(t) - 1):
        y[i+1] = rk4_step(t[i], y[i], dt, case)
        # 增加检查点，如果数值发散则提前停止
        if np.any(np.isnan(y[i+1])):
            print(f"    Calculation stopped early at t={t[i+1]:.3f} due to invalid values.")
            # 将后续结果填充为 NaN，这样绘图时曲线会在此处中断
            y[i+1:] = np.nan
            break

    print(f"  ...Done.")
    return t, y

# --- 绘图函数 ---
def plot_n_r(case):
    """绘制给定情景下不同初始条件的中子密度变化曲线"""
    plt.figure(figsize=(10, 5))
    initial_n_r_values = [0.2, 0.5, 1.0]

    for n_r0 in initial_n_r_values:
        t, y = solve(case, n_r0)
        plt.plot(t, y[:, 0], label=f'$n_{{r0}}={n_r0}$')

    plt.xlabel('时间 (s)')
    plt.ylabel('中子相对密度 $n_r(t)$')

    titles = {
        1: r'反应性正向阶跃变化 ($\Delta\rho = +0.001$)',
        2: r'反应性负向阶跃变化 ($\Delta\rho = -0.001$)',
        3: r'反应性正弦变化 ($\rho(t) = -0.001 + 0.0002 \sin(t-10)$)'
    }
    plt.title(titles.get(case, f'情景 {case}'))
    plt.legend()
    plt.grid(True)
    # 根据数值范围自动调整 Y 轴可能更合适，如果数值增长很快，手动设置 ylim 可能不佳
    # plt.ylim(bottom=0) # 可以保留确保从0开始，但上限自动
    plt.yscale('linear') # 确保是线性刻度，如果数值范围很大也可尝试 'log'
    plt.tight_layout()
    plt.show()

# --- 主程序 ---
# 运行并绘制三种情景的结果
for case in [1, 2, 3]:
    plot_n_r(case)