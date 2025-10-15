# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os

# 配置路径
ALN_RESULTS = "D:/tools/data/GII_GIX_foldseek/aln_results_rmsd.tsv"
OUTPUT_DIR = "D:/tools/data/GII_GIX_foldseek/results/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Foldseek输出列顺序
COLUMN_NAMES = ['query', 'target', 'qaln', 'taln', 'evalue', 'rmsd']

# 读取Foldseek比对结果
print("Loading Foldseek alignment results...")
try:
    aln_df = pd.read_csv(ALN_RESULTS, sep='\t', header=None, names=COLUMN_NAMES)
    print(f"Total alignments loaded: {len(aln_df)}")
except Exception as e:
    print(f"Error loading file: {e}")
    exit()

# 提取所有GII_pdb和GIX_pdb的ID
GII_ids = set()
GIX_ids = set()

for _, row in aln_df.iterrows():
    if row['query'].startswith('GII_pdb'):
        GII_ids.add(row['query'])
    if row['target'].startswith('GII_pdb'):
        GII_ids.add(row['target'])
    if row['query'].startswith('GIX_pdb'):
        GIX_ids.add(row['query'])
    if row['target'].startswith('GIX_pdb'):
        GIX_ids.add(row['target'])

print(f"Found {len(GII_ids)} unique GII_pdb structures")
print(f"Found {len(GIX_ids)} unique GIX_pdb structures")
print(f"Theoretical maximum pairs: {len(GII_ids) * len(GIX_ids)}")

# 筛选GII_pdb vs GIX_pdb的比对对
print("\nFiltering GII_pdb vs GIX_pdb alignments...")

filtered_data = []
found_pairs = set()

for idx, row in aln_df.iterrows():
    query = row['query']
    target = row['target']
    
    # 检查是否满足条件：一个以GII_pdb开头，另一个以GIX_pdb开头
    if (query.startswith('GII_pdb') and target.startswith('GIX_pdb')) or \
       (query.startswith('GIX_pdb') and target.startswith('GII_pdb')):
        # 创建标准化的键（按字母顺序排序，避免重复）
        pair_key = tuple(sorted([query, target]))
        if pair_key not in found_pairs:
            found_pairs.add(pair_key)
            filtered_data.append(row)

# 创建过滤后的DataFrame
filtered_df = pd.DataFrame(filtered_data)

print(f"Found {len(filtered_df)} unique GII_pdb vs GIX_pdb alignments")
print(f"Coverage: {len(filtered_df)/ (len(GII_ids) * len(GIX_ids)) * 100:.2f}% of theoretical maximum")

# 如果没有找到足够的比对对，显示一些示例数据来帮助诊断
if len(filtered_df) < len(GII_ids) * len(GIX_ids):
    print("\nNot all possible pairs were found. This could be due to:")
    print("1. FoldSeek filtering based on e-value or other criteria")
    print("2. Poor quality structures that couldn't be aligned")
    print("3. Memory or computational limitations during FoldSeek run")
    
    # 保存诊断信息
    with open(os.path.join(OUTPUT_DIR, 'coverage_analysis.txt'), 'w', encoding='utf-8') as f:
        f.write("Coverage Analysis\n")
        f.write("=================\n")
        f.write(f"Unique GII_pdb structures: {len(GII_ids)}\n")
        f.write(f"Unique GIX_pdb structures: {len(GIX_ids)}\n")
        f.write(f"Theoretical maximum pairs: {len(GII_ids) * len(GIX_ids)}\n")
        f.write(f"Actual pairs found: {len(filtered_df)}\n")
        f.write(f"Coverage: {len(filtered_df)/ (len(GII_ids) * len(GIX_ids)) * 100:.2f}%\n")
        f.write("\nGII_pdb structures:\n")
        for i, id in enumerate(sorted(GII_ids)):
            f.write(f"{i+1}. {id}\n")
        f.write("\nGIX_pdb structures:\n")
        for i, id in enumerate(sorted(GIX_ids)):
            f.write(f"{i+1}. {id}\n")
    
    print("Coverage analysis saved to coverage_analysis.txt")

# 提取RMSD值并转换为float
rmsd_values = filtered_df['rmsd'].astype(float).values

# 计算总体统计量
mean_rmsd = np.mean(rmsd_values)
std_rmsd = np.std(rmsd_values)
median_rmsd = np.median(rmsd_values)
min_rmsd = np.min(rmsd_values)
max_rmsd = np.max(rmsd_values)

# 输出标准RMSD统计结果
print("\n=== STANDARD RMSD STATISTICS ===")
print(f"Number of alignments: {len(rmsd_values)}")
print(f"Mean RMSD: {mean_rmsd:.4f} Å")
print(f"Std RMSD: {std_rmsd:.4f} Å")
print(f"Median RMSD: {median_rmsd:.4f} Å")
print(f"Min RMSD: {min_rmsd:.4f} Å")
print(f"Max RMSD: {max_rmsd:.4f} Å")
print("================================")

# 保存统计结果到文件
stats_file = os.path.join(OUTPUT_DIR, 'GII_vs_GIX_standard_rmsd_stats.txt')
with open(stats_file, 'w', encoding='utf-8') as f:
    f.write("Standard RMSD Statistics for GII_pdb vs GIX_pdb comparisons\n")
    f.write("============================================================\n")
    f.write(f"Number of alignments: {len(rmsd_values)}\n")
    f.write(f"Mean RMSD: {mean_rmsd:.4f} Å\n")
    f.write(f"Standard deviation RMSD: {std_rmsd:.4f} Å\n")
    f.write(f"Median RMSD: {median_rmsd:.4f} Å\n")
    f.write(f"Minimum RMSD: {min_rmsd:.4f} Å\n")
    f.write(f"Maximum RMSD: {max_rmsd:.4f} Å\n")
    f.write("\nPairwise comparisons:\n")
    f.write("Query,Target,RMSD\n")
    for _, row in filtered_df.iterrows():
        f.write(f"{row['query']},{row['target']},{row['rmsd']}\n")

print(f"Saved standard RMSD statistics to {stats_file}")

print("\nAnalysis completed successfully!")