import numpy as np
from scipy.interpolate import interp1d
from scipy.linalg import lu
import matplotlib.pyplot as plt
import matplotlib as mpl

# 设置 Matplotlib 支持中文显示
mpl.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
mpl.rcParams['axes.unicode_minus'] = False  # 解决保存图像时负号'-'显示为方块的问题

# --- 数据定义 ---
power_levels = np.array([3.2, 4.1, 9.5, 24.2, 30, 50, 100])
Tn_data = np.array([5.14, 8.00, 9.00, 6.29, 5.71, 5.71, 5.71])
Fg_data = np.array([13.00, 18.00, 10.00, 4.00, 4.00, 4.00, 4.00])
Th_data = np.array([24.29, 8.00, 4.29, 1.43, 1.14, 0.71, 0.71])
tau_data = np.array([1.43, 1.43, 1.43, 4.29, 4.29, 4.29, 4.29])
Tg = 1.429  # 常数

params_data = {
    'Tn': Tn_data,
    'Fg': Fg_data,
    'Th': Th_data,
    'tau': tau_data
}
param_names = ['Tn', 'Fg', 'Th', 'tau']

print("--- 二、(60 分)蒸汽发生器水位控制系统问题求解 ---")

# --- 问题 1: 最小二乘多项式拟合 ---
print("\n--- 问题 1: 最小二乘多项式拟合 ---")
power_eval_fit = np.array([10.0, 40.0, 80.0])
fitted_params = {}
fitting_errors = {}

plt.figure(figsize=(8, 6))
plot_index = 1

for name, data in params_data.items():
    # 构造三次多项式拟合
    coeffs = np.polyfit(power_levels, data, 3)
    poly_func = np.poly1d(coeffs)
    fitted_params[name] = poly_func

    # 计算指定功率下的近似值
    approx_values = poly_func(power_eval_fit)
    print(f"\n参数 {name} 的三次拟合多项式:")
    print(poly_func)
    print(f"在功率 {power_eval_fit}% 时的近似值: {approx_values}")

    # 计算拟合误差 (在原始数据点上的误差)
    fitted_values_at_data_points = poly_func(power_levels)
    errors = data - fitted_values_at_data_points
    sse = np.sum(errors**2) # 残差平方和
    fitting_errors[name] = {'residuals': errors, 'sse': sse}
    print(f"在原始数据点上的拟合残差: {errors}")
    print(f"残差平方和 (SSE): {sse:.4f}")

    # 绘制拟合曲线
    plt.subplot(2, 2, plot_index)
    p_smooth = np.linspace(power_levels.min(), power_levels.max(), 200)
    plt.plot(p_smooth, poly_func(p_smooth), label=f'{name} 三次拟合曲线')
    plt.scatter(power_levels, data, color='red', label=f'{name} 原始数据点')
    plt.scatter(power_eval_fit, approx_values, color='green', marker='x', s=100, label=f'{name} 拟合预测点 ({power_eval_fit}%)')
    plt.title(f'{name} 参数随功率变化的拟合曲线')
    plt.xlabel('功率水平 (%)')
    plt.ylabel(f'{name} 值')
    plt.legend()
    plt.grid(True)
    plot_index += 1

plt.tight_layout()
plt.show()


# --- 问题 2: 插值计算 ---
print("\n--- 问题 2: 插值计算 ---")
power_eval_interp = np.array([20.0, 60.0, 90.0])
interp_results_linear = {}
interp_results_spline = {}

print(f"\n使用分段线性插值和分段三次样条插值计算功率为 {power_eval_interp}% 时的参数近似值:")

for name, data in params_data.items():
    # 分段线性插值
    # 对于插值，如果插值点超出原始数据范围，需要处理边界情况
    # fill_value="extrapolate" 允许外插
    linear_interp_func = interp1d(power_levels, data, kind='linear', fill_value="extrapolate")
    linear_values = linear_interp_func(power_eval_interp)
    interp_results_linear[name] = linear_values
    print(f"\n参数 {name}:")
    print(f"  分段线性插值结果: {linear_values}")

    # 分段三次样条插值
    # 'cubic' 指的是三次样条插值
    spline_interp_func = interp1d(power_levels, data, kind='cubic', fill_value="extrapolate")
    spline_values = spline_interp_func(power_eval_interp)
    interp_results_spline[name] = spline_values
    print(f"  分段样条插值结果: {spline_values}")

# --- 问题 3: LU 分解 ---
print("\n--- 问题 3: LU 分解 ---")
power_eval_lu = np.array([30.0, 50.0, 100.0])

# 从原始数据中直接获取这些功率水平下的参数值
lu_params = {}
for p_val in power_eval_lu:
    idx = np.where(power_levels == p_val)[0][0] # 找到对应功率的索引
    lu_params[p_val] = {
        'Tn': Tn_data[idx],
        'Th': Th_data[idx],
        'tau': tau_data[idx]
        # Fg 不需要，因为它不出现在 A 矩阵中
        # Tg 是常数
    }

def get_A_matrix(Tn, Th, Tg, tau):
    """根据参数计算矩阵 A"""
    if Tn == 0 or Th == 0 or Tg == 0 or tau == 0:
        print(f"警告: 参数为零 (Tn={Tn}, Th={Th}, Tg={Tg}, tau={tau})，可能导致除零错误。")
        # 可以选择返回错误、None 或一个特殊矩阵
        return None # 或者根据情况处理

    return np.array([
        [0, 0, 0, 1/Tn if Tn != 0 else np.inf],
        [0, -1/Th if Th != 0 else np.inf, 0, -1/Tn if Tn != 0 else np.inf],
        [0, 0, -1/Tg if Tg != 0 else np.inf, 0],
        [0, 0, 0, -1/tau if tau != 0 else np.inf]
    ])

print("\n计算 A 矩阵的 LU 分解:")
for p_val in power_eval_lu:
    params_at_p = lu_params[p_val]
    Tn_val = params_at_p['Tn']
    Th_val = params_at_p['Th']
    tau_val = params_at_p['tau']

    print(f"\n功率水平 = {p_val}%:")
    print(f"  参数: Tn={Tn_val}, Th={Th_val}, Tg={Tg}, tau={tau_val}")

    A = get_A_matrix(Tn_val, Th_val, Tg, tau_val)

    if A is not None:
        print("  矩阵 A:")
        print(A)

        # 执行 LU 分解 (P @ A = L @ U)
        # P 是置换矩阵, L 是下三角矩阵, U 是上三角矩阵
        try:
            P, L, U = lu(A)
            print("  LU 分解结果 (P @ A = L @ U):")
            print("  置换矩阵 P:")
            print(P)
            print("  下三角矩阵 L:")
            print(L)
            print("  上三角矩阵 U:")
            print(U)

            # 验证分解 P @ A 约等于 L @ U
            # print("  验证 P @ A:")
            # print(P @ A)
            # print("  验证 L @ U:")
            # print(L @ U)
            # print(f"  分解是否接近: {np.allclose(P @ A, L @ U)}")

        except np.linalg.LinAlgError as e:
             print(f"  LU 分解失败: {e}") # 例如，如果矩阵是奇异的
        except ZeroDivisionError:
             print("  LU 分解失败: 遇到除零错误（参数可能为零）")
    else:
        print("  无法计算矩阵 A (参数为零?)")

# 验证 LU 分解是否正确
for p_val in power_eval_lu:
    params_at_p = lu_params[p_val]
    Tn_val = params_at_p['Tn']
    Th_val = params_at_p['Th']
    tau_val = params_at_p['tau']

    print(f"\n功率水平 = {p_val}%:")
    print(f"  参数: Tn={Tn_val}, Th={Th_val}, Tg={Tg}, tau={tau_val}")

    A = get_A_matrix(Tn_val, Th_val, Tg, tau_val)

    if A is not None:
        print("  矩阵 A:")
        print(A)

        # 执行 LU 分解 (P @ A = L @ U)
        try:
            P, L, U = lu(A)
            print("  LU 分解结果 (P @ A = L @ U):")
            print("  置换矩阵 P:")
            print(P)
            print("  下三角矩阵 L:")
            print(L)
            print("  上三角矩阵 U:")
            print(U)

            # 验证 P @ A 是否接近 L @ U
            PA = np.dot(P, A)
            LU = np.dot(L, U)
            print(f"  P @ A: \n{PA}")
            print(f"  L @ U: \n{LU}")
            print(f"  P @ A 是否接近 L @ U: {np.allclose(PA, LU)}")

        except np.linalg.LinAlgError as e:
             print(f"  LU 分解失败: {e}") # 例如，如果矩阵是奇异的
        except ZeroDivisionError:
             print("  LU 分解失败: 遇到除零错误（参数可能为零）")
    else:
        print("  无法计算矩阵 A (参数为零?)")

print("\n--- 计算完成 ---")