import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# 设置全局字体为支持中文的字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 使用 SimHei 字体（黑体）
matplotlib.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 参数
beta = 0.006019
Lambda = 0.00002
lambda_ = 0.150
lambda_I = 2.9e-5
lambda_X = 2.1e-5
theta = 0.3898
vartheta = 2.807

# 反应性函数
def rho_t(t, case):
    if case == 1:  # 正向阶跃
        return 0 if t < 10 else 0.001
    elif case == 2:  # 负向阶跃
        return 0 if t < 10 else -0.001
    elif case == 3:  # 正弦变化
        return 0 if t < 10 else -0.001 + 0.0002 * np.sin(t - 10)

# 微分方程
def derivatives(t, y, case):
    n_r, c_r, I_r, X_r = y
    rho = rho_t(t, case)
    dn_r_dt = ((rho - beta) / Lambda) * n_r + (beta / Lambda) * c_r
    dc_r_dt = lambda_ * (n_r - c_r)
    dI_r_dt = lambda_I * (n_r - I_r)
    term1 = (I_r + theta * n_r) / (1 + theta)
    term2 = (1 + vartheta * n_r) / (1 + theta) * X_r
    dX_r_dt = lambda_X * (1 + vartheta) * (term1 - term2)
    return [dn_r_dt, dc_r_dt, dI_r_dt, dX_r_dt]

# 四阶龙格-库塔步进
def rk4_step(t, y, dt, case):
    k1 = np.array(derivatives(t, y, case))
    k2 = np.array(derivatives(t + dt/2, y + dt/2 * k1, case))
    k3 = np.array(derivatives(t + dt/2, y + dt/2 * k2, case))
    k4 = np.array(derivatives(t + dt, y + dt * k3, case))
    return y + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)

# 求解函数
def solve(case, n_r0, t_max=400, dt=0.1):
    t = np.arange(0, t_max + dt, dt)
    y = np.zeros((len(t), 4))
    y[0] = [n_r0, n_r0, n_r0, n_r0 * (1 + theta) / (1 + vartheta * n_r0)]
    for i in range(1, len(t)):
        y[i] = rk4_step(t[i-1], y[i-1], dt, case)
    return t, y

# 绘制函数
def plot_n_r(case):
    plt.figure(figsize=(10, 6))
    for n_r0 in [0.2, 0.5, 1.0]:
        t, y = solve(case, n_r0)
        plt.plot(t, y[:, 0], label=f'n_r0={n_r0}')
    plt.xlabel('时间 (s)')
    plt.ylabel('中子相对密度 n_r(t)')
    titles = {1: '反应性正向阶跃变化', 2: '反应性负向阶跃变化', 3: '反应性正弦变化'}
    plt.title(titles[case])
    plt.legend()
    plt.grid(True)
    plt.show()

# 主程序
for case in [1, 2, 3]:
    plot_n_r(case)