import os
import glob
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
import matplotlib.pyplot as plt
import seaborn as sns

# 设置路径
data_dir = "path.."
aa_fasta_path = "path.."
pdb_files = glob.glob(os.path.join(data_dir, "*.pdb"))

# 读取对齐的氨基酸序列（用来确定序列长度和位置）
from Bio import SeqIO

# 读取所有记录
records = list(SeqIO.parse(aa_fasta_path, "fasta"))

# 如果你只关心第一条
aa_record = records[0]
aa_sequence = str(aa_record.seq)
seq_length = len(aa_sequence)

# 初始化存储pLDDT值的矩阵（结构数×序列长度）
plddt_values = []

# 初始化PDB解析器
parser = PDBParser(QUIET=True)

for pdb_file in pdb_files:
    structure_id = os.path.basename(pdb_file).split('.')[0]
    structure = parser.get_structure(structure_id, pdb_file)
    # 假设结构中只有一条链，或者目标链为链A
    model = structure[0]
    chain = next(model.get_chains())

    # 提取每个残基的pLDDT (B因子)
    # 按照对齐位置提取
    residues = list(chain.get_residues())
    res_bfactors = []
    for res in residues:
        # 只考虑标准残基
        if res.id[0] == ' ':
            res_bfactors.append(res['CA'].bfactor)
    # 填补到对齐长度，缺失位置用np.nan
    res_bfactors_full = np.full(seq_length, np.nan)
    for i, res in enumerate(residues):
        if i >= seq_length:
            break
        res_bfactors_full[i] = res['CA'].bfactor
    plddt_values.append(res_bfactors_full)

# 将数据转换为DataFrame，行是结构，列是位置
df_plddt = pd.DataFrame(plddt_values)

# 计算每个位置的平均值（跳过NaN）
mean_plddt = df_plddt.mean(axis=0)

# 绘制热图
plt.figure(figsize=(20, 4))
sns.heatmap([mean_plddt], cmap='rainbow', cbar=True, vmin=0, vmax=100, xticklabels=100)
plt.xlabel("Alignment Position")
plt.title("Average pLDDT per Position")
plt.show()