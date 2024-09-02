import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# 读取包含SMILES的CSV文件
df = pd.read_csv('smile.csv')

# 创建一个空列表，用于存储分子对象
molecules = []

# 遍历CSV中的每一行，将SMILES字符串转换为分子对象并添加到列表中
for smile in df['SMILES']:
    molecule = Chem.MolFromSmiles(smile)
    if molecule is not None:
        molecules.append(molecule)

# 设置绘图参数，可以根据需要进行调整
draw_options = Draw.DrawingOptions()
draw_options.atomLabelFontSize = 30
draw_options.dotsPerAngstrom = 100
draw_options.bondLineWidth = 3.0

# 绘制每个分子的分子图，并保存为图片文件
for idx, mol in enumerate(molecules):
    Draw.MolToFile(mol, f'molecule_{idx}.png', size=(400, 400), options=draw_options)
