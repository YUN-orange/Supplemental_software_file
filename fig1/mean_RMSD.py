# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 加载数据
df = pd.read_csv('D:/tools/data/GII.3_foldseek/results/residue_rmsd_contributions.csv')

# 数据预处理
df = df.dropna(subset=['rmsd_contribution'])
conservation = df.groupby('residue_number')['rmsd_contribution'].agg(['mean', 'std'])

# 创建图形并设置全局字体
plt.rcParams.update({'font.size': 16})  # 设置全局基础字体大小
plt.figure(figsize=(17, 8), dpi=300)

# 绘制折线图 (均值) - 移除label参数
plt.plot(
    conservation.index,
    conservation['mean'],
    color='royalblue',
    linewidth=2.5  # 增加线宽
)

# 绘制标准差阴影区域 - 移除label参数
plt.fill_between(
    x=conservation.index,
    y1=conservation['mean'] - conservation['std'],
    y2=conservation['mean'] + conservation['std'],
    color='skyblue',
    alpha=0.4
)

# 添加全局均值参考线 (保留label)
global_mean = df['rmsd_contribution'].mean()
plt.axhline(
    y=global_mean,
    color='red',
    linestyle='--',
    linewidth=2.0,  # 增加线宽
    alpha=0.8,
    label=f'Global Mean ({global_mean:.2f}Å)'
)

# 设置图形属性
plt.title('Residue RMSD Stability Analysis', pad=20, fontsize=40)
plt.xlabel('Residue Number', fontsize=38)
plt.ylabel('RMSD Contribution (Å)', fontsize=38)
plt.grid(alpha=0.2, linestyle=':')

# 自动设置纵坐标范围和刻度（保留一位小数）
max_value = conservation['mean'].max() + conservation['std'].max()
upper_limit = np.ceil(max_value * 10) / 10  # 向上取整到一位小数
yticks = np.arange(0, upper_limit + 0.1, 0.5)  # 每0.5一个刻度
plt.ylim(0, upper_limit)
plt.yticks(yticks, fontsize=36)

# 设置横坐标刻度字体大小
plt.xticks(fontsize=36)

# 智能调整X轴刻度
if len(conservation) > 50:
    step = max(1, len(conservation) // 20)
    plt.xticks(conservation.index[::step], rotation=45)
else:
    plt.xticks(conservation.index, rotation=45)

# 添加图例 (现在只会显示Global Mean)
plt.legend(
    loc='best', 
    frameon=True,
    fontsize=36,
    title_fontsize=38
)

# 优化布局并保存
plt.tight_layout()
plt.savefig(path.., 
            bbox_inches='tight', transparent=False)
plt.close()

print("Simplified RMSD stability plot with only Global Mean legend saved successfully.")