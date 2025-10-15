import os
import numpy as np
import subprocess

# 配置路径
pdb_filename = "AB039776.pdb"  # 你的参考PDB文件名
pdb_dir = "D:/tools/data/GII_pdbs"
pdb_path = os.path.join(pdb_dir, pdb_filename)

# 你的对齐序列长度（请确认）
# 假设seq_length就是你之前定义的
seq_length = len([323])  # 替换为你的序列长度

# 你的平均pLDDT值（请确保在之前的代码中已计算）
# 示例：只演示用随机数作为占位
# 在实际应用中，请用你之前计算的mean_plddt值
mean_plddt = np.random.uniform(0, 100, size=seq_length)  # 替换为你的真实数据

# 生成PyMOL脚本路径
pymol_path = r"C:\Program Files\Pymol\PyMOLWin.exe"
pml_path = "D:/tools/data/GII_pdbs/mapping_plddt.pml"

# 生成PyMOL脚本内容
with open(pml_path, "w") as f:
    f.write(f'load {pdb_path}, protein\n')
    # 设置每个残基的b因子为对应的pLDDT
    for resi in range(1, seq_length + 1):
        plddt_value = mean_plddt[resi - 1]
        # 忽略NaN或无效值
        if not np.isnan(plddt_value):
            f.write(f'alter (resi {resi}), b={plddt_value}\n')
    f.write('sort\n')
    f.write('spectrum b, rainbow, all, minimum=0, maximum=100\n')
    f.write('set ray_shadows, 0\n')
    f.write('set ray_opaque_background, 0\n')
    f.write('bg_color white\n')
    # 保存图片路径
    output_image = os.path.join(pdb_dir, "structure_plddt_mapped.png")
    f.write(f'png {output_image}, width=800, height=600, dpi=300\n')
    f.write('quit\n')

# 调用PyMOL执行脚本
subprocess.run([pymol_path, "-c", pml_path])

print(f"映射完成，图片已保存到：{output_image}")