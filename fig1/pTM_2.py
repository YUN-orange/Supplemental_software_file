import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 设置中文显示（如果需要）
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False    # 用来正常显示负号

# 设置输出目录
output_dir = r"D:\tools\data\ptm_analysis\visualization"
os.makedirs(output_dir, exist_ok=True)

# 读取pTM数据
ptm_df = pd.read_csv(r"D:\tools\data\ptm_analysis\ptm_summary.csv")

# 基本分布可视化
def plot_basic_distribution(df):
    """绘制pTM值的基本分布图"""
    plt.figure(figsize=(15, 8))
    
    # 直方图
    plt.subplot(1, 2, 1)
    sns.histplot(df['ptm'], bins=30, kde=True, color='royalblue')
    plt.axvline(x=0.7, color='r', linestyle='--', label='高置信度 (0.7)')
    plt.axvline(x=0.5, color='orange', linestyle='--', label='中等置信度 (0.5)')
    plt.xlabel('pTM值')
    plt.ylabel('数量')
    plt.title('pTM值分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 箱线图
    plt.subplot(1, 2, 2)
    sns.boxplot(y=df['ptm'], color='lightgreen')
    plt.ylabel('pTM值')
    plt.title('pTM值箱线图')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ptm_basic_distribution.png'), dpi=300)
    plt.close()
    print("基本分布图已保存")

# 执行可视化
print("\n开始可视化分析...")
plot_basic_distribution(ptm_df)
print("\n所有可视化分析已完成!")
print(f"结果保存在: {output_dir}")