import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# 设置全局字体为支持中文的字体
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 使用 SimHei 字体（黑体）
matplotlib.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# -------------------
# 参数设置
# -------------------
beta = 0.006019
Lambda = 0.00002
lam = 0.150
lam_I = 2.9e-5
lam_X = 2.1e-5
gamma_I = 0.059
gamma_X = 0.023
sigma_X = 3.5e-18
P0 = 3000  # MW
E_ff = 3.2e-11  # MW·s
D = 316  # cm
h = 355  # cm

V = 2 * np.pi * (0.5 * D)**2 * h  # 体积，单位：cm³
phi0 = P0 / (E_ff * V)  # 中子注量率
theta = gamma_X / gamma_I
vartheta = sigma_X * phi0 / lam_X

# -------------------
# 初值函数
# -------------------
def initial_conditions(n_r0):
    c_r0 = n_r0
    I_r0 = n_r0
    X_r0 = (I_r0 + theta * n_r0) / (1 + vartheta * n_r0)
    return np.array([n_r0, c_r0, I_r0, X_r0])

# -------------------
# 反应性函数
# -------------------
def rho_step_pos(t):
    return 0 if t < 10 else 0.001

def rho_step_neg(t):
    return 0 if t < 10 else -0.001

def rho_sin(t):
    return 0 if t < 10 else -0.001 + 0.0002 * np.sin(t - 10)

# -------------------
# 微分方程组定义
# -------------------
def derivatives(t, y, rho_func):
    n_r, c_r, I_r, X_r = y
    rho = rho_func(t)
    dn_dt = ((rho - beta) / Lambda) * n_r + (beta / Lambda) * c_r
    dc_dt = lam * (n_r - c_r)
    dI_dt = lam_I * (n_r - I_r)
    dX_dt = lam_X * (1 + vartheta) * ((I_r + theta * n_r) / (1 + theta) - ((1 + vartheta * n_r) / (1 + theta)) * X_r)
    return np.array([dn_dt, dc_dt, dI_dt, dX_dt])

# -------------------
# RK4积分器
# -------------------
def runge_kutta4(f, y0, t_span, dt, rho_func):
    t0, t_end = t_span
    ts = np.arange(t0, t_end + dt, dt)
    ys = np.zeros((len(ts), len(y0)))
    ys[0] = y0
    for i in range(1, len(ts)):
        t = ts[i - 1]
        y = ys[i - 1]
        k1 = f(t, y, rho_func)
        k2 = f(t + dt/2, y + dt/2 * k1, rho_func)
        k3 = f(t + dt/2, y + dt/2 * k2, rho_func)
        k4 = f(t + dt, y + dt * k3, rho_func)
        ys[i] = y + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    return ts, ys

# -------------------
# 作图函数
# -------------------
def plot_result(ts, ys, label):
    plt.plot(ts, ys[:, 0], label=label)  # 只画 n_r

# -------------------
# 主程序：三种初值 & 三种rho情况
# -------------------
dt = 0.1
t_span = (0, 400)
rho_funcs = [rho_step_pos, rho_step_neg, rho_sin]
rho_labels = ["正阶跃", "负阶跃", "正弦变化"]
n_r0_list = [0.2, 0.5, 1.0]

for i, rho_func in enumerate(rho_funcs):
    plt.figure(figsize=(8, 5))
    for n_r0 in n_r0_list:
        y0 = initial_conditions(n_r0)
        ts, ys = runge_kutta4(derivatives, y0, t_span, dt, rho_func)
        plot_result(ts, ys, f"n_r0={n_r0}")
    plt.title(f"中子相对密度变化 - {rho_labels[i]}")
    plt.xlabel("时间 t (s)")
    plt.ylabel("n_r(t)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()