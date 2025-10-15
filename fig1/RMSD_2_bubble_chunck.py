# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

# 初始化用于累积统计量的数据结构
sums = {}
squares = {}
counts = {}

# 分块读取和处理数据
chunk_size = 1000000  # 根据内存调整块大小
chunk_iterator = pd.read_csv(
    'D:/tools/data/GII.4_foldseek/rmsd_results/residue_rmsd_contributions.csv',
    chunksize=chunk_size,
    usecols=['residue_number', 'rmsd_contribution']  # 只读取需要的列
)

for i, chunk in enumerate(chunk_iterator):
    print(f"Processing chunk {i+1}")
    
    # 删除缺失值
    chunk = chunk.dropna(subset=['rmsd_contribution'])
    
    # 计算每个残基的统计量
    chunk_grouped = chunk.groupby('residue_number')['rmsd_contribution'].agg(['sum', 'count'])
    chunk_grouped['sum_of_squares'] = chunk.groupby('residue_number')['rmsd_contribution'].apply(lambda x: (x**2).sum())
    
    # 累积统计量
    for residue, row in chunk_grouped.iterrows():
        if residue not in sums:
            sums[residue] = 0
            squares[residue] = 0
            counts[residue] = 0
        
        sums[residue] += row['sum']
        squares[residue] += row['sum_of_squares']
        counts[residue] += row['count']

# 计算最终的均值和标准差
conservation_data = []
for residue in sums:
    mean = sums[residue] / counts[residue]
    variance = (squares[residue] / counts[residue]) - (mean ** 2)
    std = np.sqrt(variance) if variance > 0 else 0
    cv = std / mean if mean != 0 else 0
    
    conservation_data.append({
        'residue_number': residue,
        'mean': mean,
        'std': std,
        'cv': cv
    })

# 创建DataFrame
conservation = pd.DataFrame(conservation_data).set_index('residue_number')

# 计算全局均值
global_mean = sum(sums.values()) / sum(counts.values())

# 创建图形并设置全局字体大小
plt.rcParams.update({
    'font.size': 18,
    'legend.fontsize': 16
})
plt.figure(figsize=(17, 10), dpi=300)  # 增加高度到10英寸

# 设置颜色映射和标准化
norm = Normalize(vmin=conservation['cv'].min(), vmax=conservation['cv'].max())
cmap = cm.get_cmap('coolwarm')
sm = cm.ScalarMappable(norm=norm, cmap=cmap)

# 动态计算气泡大小
min_size = 20
max_size = 300
size_range = conservation['std'].max() - conservation['std'].min()
if size_range > 0:
    bubble_sizes = min_size + (conservation['std'] - conservation['std'].min()) / size_range * (max_size - min_size)
else:
    bubble_sizes = np.full_like(conservation['std'], (min_size + max_size)/2)

# 绘制散点图
scatter = plt.scatter(
    x=conservation.index,
    y=conservation['mean'],
    c=conservation['cv'],
    s=bubble_sizes,
    cmap='coolwarm',
    alpha=0.7,
    edgecolors='black',
    linewidths=0.5
)

# 添加参考线
plt.axhline(
    y=global_mean,
    color='red',
    linestyle='--',
    linewidth=2,
    alpha=0.7,
    label=f'Global Mean ({global_mean:.2f}Å)'
)

# 添加颜色条
cbar = plt.colorbar(scatter, pad=0.01)
cbar.set_label('Coefficient of Variation (CV)', rotation=270, labelpad=25, fontsize=18)
cbar.ax.tick_params(labelsize=18)

# 添加气泡大小图例
for size in [min_size, (min_size + max_size)//2, max_size]:
    # 计算对应的标准差范围
    if size_range > 0:
        std_value = conservation['std'].min() + (size - min_size) / (max_size - min_size) * size_range
    else:
        std_value = conservation['std'].mean()
    
    plt.scatter(
        [], [],
        s=size,
        c='gray',
        alpha=0.5,
        edgecolors='black',
        label=f'Std={std_value:.2f}Å'
    )

# 设置图形属性
plt.title(
    'Residue Conservation Analysis\n'
    'Color: Coefficient of Variation | Bubble Size: Standard Deviation',
    pad=20,
    fontsize=18,
    y=1.02  # 稍微上移标题
)
plt.xlabel('Residue Number', fontsize=18)
plt.ylabel('Mean RMSD Contribution (Å)', fontsize=18)
plt.grid(alpha=0.2, linestyle=':')

# 设置坐标轴刻度
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# 智能调整X轴刻度
if len(conservation) > 50:
    step = max(1, len(conservation) // 20)
    plt.xticks(conservation.index[::step], rotation=45)
else:
    plt.xticks(conservation.index, rotation=45)

# 添加并调整图例位置
legend = plt.legend(
    loc='upper center',
    bbox_to_anchor=(0.5, 0.98),  # 向下调整位置
    ncol=4,
    frameon=True,
    handletextpad=1.5,
    columnspacing=2.0,
    markerscale=0.8,
    borderaxespad=0.5
)

# 调整布局
plt.subplots_adjust(top=0.88)  # 调整顶部空间

# 保存图形
plt.tight_layout()
plt.savefig(
    'D:/tools/data/GII.4_foldseek/rmsd_results/conservation_analysis_final.png',
    bbox_inches='tight',
    transparent=False,
    dpi=300
)
plt.close()

print("Final conservation analysis plot with optimized legend position saved successfully.")