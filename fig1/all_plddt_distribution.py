import os
from Bio.PDB import PDBParser
import pandas as pd
import numpy as np

# 设置路径
input_dir = r"D:\tools\data\input_pdbs"
output_dir = r"D:\tools\data\plddt_analysis"  # 专门定义输出目录

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 输出文件路径
output_csv = os.path.join(output_dir, "all_plddt_data.csv")
summary_csv = os.path.join(output_dir, "structure_summary.csv")

all_data = []
summary_list = []
error_log = []

parser = PDBParser()
processed_files = 0

print(f"开始处理目录: {input_dir} 中的PDB文件...")
print(f"共找到 {len(os.listdir(input_dir))} 个文件")

for pdb_file in os.listdir(input_dir):
    if pdb_file.endswith(".pdb"):
        filepath = os.path.join(input_dir, pdb_file)
        try:
            structure = parser.get_structure(pdb_file[:-4], filepath)
            model = structure[0]
            
            print(f"处理文件中: {pdb_file}...")
            
            for chain in model:
                chain_id = chain.id
                chain_data = []
                for residue in chain:
                    # 只处理标准氨基酸残基
                    if residue.id[0] == ' ':
                        res_num = residue.id[1]
                        res_name = residue.get_resname()
                        
                        # 获取第一个原子的pLDDT值
                        for atom in residue:
                            plddt = atom.get_bfactor()
                            break
                        
                        chain_data.append({
                            'filename': pdb_file,
                            'chain': chain_id,
                            'residue_number': res_num,
                            'residue_name': res_name,
                            'plddt': plddt
                        })
                
                if chain_data:
                    chain_df = pd.DataFrame(chain_data)
                    all_data.append(chain_df)
                    
                    # 计算统计量
                    plddt_values = chain_df['plddt']
                    summary_list.append({
                        'filename': pdb_file,
                        'chain': chain_id,
                        'mean_plddt': np.mean(plddt_values),
                        'median_plddt': np.median(plddt_values),
                        'min_plddt': np.min(plddt_values),
                        'max_plddt': np.max(plddt_values),
                        'std_plddt': np.std(plddt_values),
                        'residues_count': len(plddt_values)
                    })
            
            processed_files += 1
            if processed_files % 10 == 0:
                print(f"已处理 {processed_files} 个文件...")
                
        except Exception as e:
            error_msg = f"解析 {pdb_file} 时出错: {str(e)}"
            print(error_msg)
            error_log.append(error_msg)

# 保存数据
if all_data:
    print("合并数据并保存...")
    master_df = pd.concat(all_data, ignore_index=True)
    master_df.to_csv(output_csv, index=False)
    print(f"已保存残基级数据到: {output_csv}")
    
    summary_df = pd.DataFrame(summary_list)
    summary_df.to_csv(summary_csv, index=False)
    print(f"已保存结构摘要数据到: {summary_csv}")
else:
    print("警告: 未提取到任何数据!")

# 保存错误日志
if error_log:
    error_log_path = os.path.join(output_dir, "error_log.txt")
    with open(error_log_path, "w") as f:
        f.write("\n".join(error_log))
    print(f"发现 {len(error_log)} 个错误，已保存到: {error_log_path}")

print("数据处理完成!")

import seaborn as sns
import matplotlib.pyplot as plt

# 1. 全局平均pLDDT分布图
plt.figure(figsize=(10, 6))
sns.histplot(summary_df['mean_plddt'], bins=30, kde=True)
plt.axvline(x=90, color='g', linestyle='--', label='Very High (90)')
plt.axvline(x=70, color='b', linestyle='--', label='High (70)')
plt.axvline(x=50, color='r', linestyle='--', label='Medium (50)')
plt.xlabel('Mean pLDDT')
plt.ylabel('Number of Structures')
plt.title('Distribution of Global Mean pLDDT Scores')
plt.legend()
plt.savefig(r"D:\tools\data\plddt_analysis\mean_plddt_distribution.png")
plt.show()

# 2. 箱线图展示统计量分布 (Mean, Min, Max)
plt.figure(figsize=(10, 6))
sns.boxplot(data=summary_df[['mean_plddt', 'min_plddt', 'max_plddt']])
plt.ylabel('pLDDT Score')
plt.title('Distribution of Mean, Min, and Max pLDDT per Structure')
plt.savefig(r"D:\tools\data\plddt_analysis\plddt_summary_boxplot.png")
plt.show()