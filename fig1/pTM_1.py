import os
import json
import pandas as pd
import re
import glob

# 设置路径
input_dir = r"D:\tools\data\GII.4_pdbs"
output_dir = r"D:\tools\data\ptm_analysis"
os.makedirs(output_dir, exist_ok=True)

# 调试信息
print(f"开始处理目录: {input_dir} 中的预测结果...")
print(f"找到 {len(os.listdir(input_dir))} 个文件/目录")

ptm_data = []  # 存储所有pTM数据
error_log = []  # 错误日志
processed_count = 0
found_count = 0

# 提取基因型的辅助函数
def extract_genotype(filename):
    """从文件名中提取基因型信息"""
    # 使用正则表达式匹配常见的诺如病毒基因型命名模式
    patterns = [
        r'G[IVXL]+\.\d+[A-Za-z]*',  # GI.1, GII.4_Sydney
        r'NoV_G[IVXL]+\.\d+',       # NoV_GI.1
        r'Norovirus_[A-Za-z]+\d+',  # Norovirus_GII4
        r'[A-Z]{2}\d+_\d+',         # GII4_2012
        r'P_Domain_[A-Za-z\d]+',    # P_Domain_GI1
        r'G[IVXL]+\d+',             # GI1, GII4
        r'[A-Z]+\d+_[A-Za-z]+'      # GII4_Sydney
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename)
        if match:
            return match.group(0)
    
    # 如果都匹配不到，返回文件名的前10个字符作为标识
    return filename[:10] if len(filename) > 10 else filename

# 搜索所有可能的pTM文件
print("搜索所有可能的pTM文件...")
ptm_files = []

# 查找所有可能的JSON文件
# 搜索所有json文件
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if file.endswith('.json'):
            ptm_files.append(os.path.join(root, file))

print(f"找到 {len(ptm_files)} 个可能的pTM文件")

# 处理每个pTM文件
for ptm_file in ptm_files:
    try:
        processed_count += 1
        
        # 获取基础信息
        dir_name = os.path.basename(os.path.dirname(ptm_file))
        file_name = os.path.basename(ptm_file)
        
        # 查找关联的PDB文件
        pdb_files = glob.glob(os.path.join(os.path.dirname(ptm_file), "*.pdb"))
        pdb_file = pdb_files[0] if pdb_files else "未找到关联PDB"
        
        # 读取JSON文件
        with open(ptm_file, 'r') as f:
            data = json.load(f)
        
        # 提取pTM值 - 适应不同版本的AlphaFold输出
        ptm_value = None
        
        # 尝试在顶层键中查找
        for key in ['ptm', 'pTM', 'predicted_tm_score', 'iptm', 'plddt']:
            if key in data:
                ptm_value = data[key]
                break
        
        # 如果没找到，尝试在模型数据中查找
        if ptm_value is None:
            for model_key in ['model_1', 'model_2', 'model_3', 'model_4', 'model_5']:
                if model_key in data:
                    model_data = data[model_key]
                    for key in ['ptm', 'pTM', 'predicted_tm_score']:
                        if key in model_data:
                            ptm_value = model_data[key]
                            break
                    if ptm_value is not None:
                        break
        
        if ptm_value is not None:
            # 添加到数据列表
            ptm_data.append({
                "source_dir": dir_name,
                "pdb_file": os.path.basename(pdb_file),
                "ptm": ptm_value,
                "source_file": file_name,
                "file_path": ptm_file
            })
            found_count += 1
            print(f"找到 pTM 值: {ptm_value} ({file_name})")
        else:
            error_msg = f"{ptm_file}: pTM值未找到"
            error_log.append(error_msg)
            print(error_msg)
    
    except Exception as e:
        error_msg = f"{ptm_file}: 解析错误 - {str(e)}"
        error_log.append(error_msg)
        print(error_msg)
    
    if processed_count % 10 == 0:
        print(f"已处理 {processed_count} 个文件，找到 {found_count} 个pTM值")

# 创建DataFrame
if ptm_data:
    ptm_df = pd.DataFrame(ptm_data)
    
    # 添加基因型信息
    if 'pdb_file' in ptm_df.columns:
        ptm_df["genotype"] = ptm_df["pdb_file"].apply(extract_genotype)
        print("成功添加基因型信息")
    else:
        print("警告: DataFrame中没有'pdb_file'列")
    
    # 保存结果
    output_csv = os.path.join(output_dir, "ptm_summary.csv")
    ptm_df.to_csv(output_csv, index=False)
    print(f"已保存 {len(ptm_df)} 条pTM记录到: {output_csv}")
    
    # 打印DataFrame信息用于调试
    print("\nDataFrame信息:")
    print(f"列名: {ptm_df.columns.tolist()}")
    print(f"前几行数据:\n{ptm_df.head()}")
else:
    print("警告: 未提取到任何pTM数据")

# 保存错误日志
if error_log:
    error_log_path = os.path.join(output_dir, "ptm_errors.log")
    with open(error_log_path, "w") as f:
        f.write("\n".join(error_log))
    print(f"发现 {len(error_log)} 个错误，已保存到: {error_log_path}")

print("处理完成!")