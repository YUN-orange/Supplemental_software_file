import os
import numpy as np
import pandas as pd
from Bio.PDB import *

# 配置路径和参数
PDB_DIR = "path.."
ALN_RESULTS = "path.."
OUTPUT_DIR = "path.."
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Foldseek输出列顺序
COLUMN_NAMES = ['query', 'target', 'qaln', 'taln', 'evalue', 'rmsd']

def load_pdb_structure(pdb_id, pdb_dir):
    """加载PDB结构并返回CA原子列表"""
    try:
        parser = PDBParser(QUIET=True)
        # 尝试不同大小写组合
        for fname in [f"{pdb_id}.pdb", f"{pdb_id.lower()}.pdb", f"{pdb_id.upper()}.pdb"]:
            pdb_path = os.path.join(pdb_dir, fname)
            if os.path.exists(pdb_path):
                structure = parser.get_structure(pdb_id, pdb_path)
                ca_atoms = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if 'CA' in residue:
                                ca_atoms.append(residue['CA'])
                return ca_atoms
        raise FileNotFoundError(f"No PDB file found for {pdb_id}")
    except Exception as e:
        print(f"Error loading PDB {pdb_id}: {str(e)}")
        return None

def calculate_residue_rmsd_contributions_batch(df_batch):
    """计算每个残基的RMSD贡献（批量处理）"""
    results = []
    
    for _, row in df_batch.iterrows():
        try:
            query_id = row['query']
            target_id = row['target']
            
            # 加载结构
            q_ca = load_pdb_structure(query_id, PDB_DIR)
            t_ca = load_pdb_structure(target_id, PDB_DIR)
            
            if q_ca is None or t_ca is None:
                continue
            
            # 提取对齐信息
            qaln, taln = row['qaln'], row['taln']
            
            # 验证比对长度
            if len(qaln) != len(taln):
                continue
            
            # 提取对齐的CA原子
            q_atoms = []
            t_atoms = []
            q_pos = t_pos = 0
            
            for qa, ta in zip(qaln, taln):
                if qa != '-' and ta != '-':
                    if q_pos < len(q_ca) and t_pos < len(t_ca):
                        q_atoms.append(q_ca[q_pos])
                        t_atoms.append(t_ca[t_pos])
                    q_pos += 1
                    t_pos += 1
                elif qa != '-':
                    q_pos += 1
                elif ta != '-':
                    t_pos += 1
            
            if not q_atoms:
                continue
            
            # 计算最佳拟合RMSD
            sup = Superimposer()
            sup.set_atoms(q_atoms, t_atoms)
            rmsd = sup.rms
            
            # 计算每个残基的贡献
            residue_contrib = {}
            q_pos = 0
            for i, (qa, ta) in enumerate(zip(qaln, taln)):
                if qa != '-' and ta != '-':
                    if q_pos < len(q_ca):
                        res_num = q_ca[q_pos].get_parent().id[1]
                        diff = q_ca[q_pos].get_coord() - t_ca[q_pos].get_coord()
                        dist_sq = np.sum(diff**2)
                        residue_contrib[res_num] = dist_sq
                    q_pos += 1
            
            # 归一化贡献
            total_dist_sq = sum(residue_contrib.values())
            scale_factor = (rmsd**2 * len(residue_contrib)) / total_dist_sq if total_dist_sq > 0 else 0
            
            for res_num in residue_contrib:
                residue_contrib[res_num] = np.sqrt(residue_contrib[res_num] * scale_factor)
            
            results.append({
                'query': query_id,
                'target': target_id,
                'total_rmsd': rmsd,
                'residue_contributions': residue_contrib,
                'aligned_length': len(q_atoms)
            })
            
        except Exception as e:
            print(f"Error processing {row['query']} vs {row['target']}: {str(e)}")
    
    return results

def main():
    # 1. 读取比对结果
    print("Loading alignment results...")
    try:
        # 使用迭代器分批读取
        chunksize = 10000  # 每批处理10,000行
        output_file = os.path.join(OUTPUT_DIR, "residue_rmsd_contributions.csv")
        
        # 写入CSV文件头
        pd.DataFrame(columns=['query', 'target', 'residue_number', 'rmsd_contribution', 'total_rmsd', 'aligned_length']).to_csv(output_file, index=False)
        
        total_processed = 0
        
        for chunk in pd.read_csv(ALN_RESULTS, sep='\t', header=None, names=COLUMN_NAMES, chunksize=chunksize):
            print(f"Processing batch of {len(chunk)} alignments...")
            
            # 2. 计算残基RMSD贡献
            results = calculate_residue_rmsd_contributions_batch(chunk)
            
            if results:
                # 3. 保存计算结果（追加到文件）
                output_data = []
                for res in results:
                    for r_num, val in res['residue_contributions'].items():
                        output_data.append({
                            'query': res['query'],
                            'target': res['target'],
                            'residue_number': r_num,
                            'rmsd_contribution': val,
                            'total_rmsd': res['total_rmsd'],
                            'aligned_length': res['aligned_length']
                        })
                
                result_df = pd.DataFrame(output_data)
                result_df.to_csv(output_file, mode='a', header=False, index=False)
                total_processed += len(results)
                print(f"Processed {total_processed} alignments so far")
            
            # 手动清理内存
            del chunk, results, output_data, result_df
            
        print(f"Total processed: {total_processed} alignments")
        print(f"Results saved to: {output_file}")
        
    except Exception as e:
        print(f"Error processing alignment file: {str(e)}")
        return

if __name__ == "__main__":
    main()