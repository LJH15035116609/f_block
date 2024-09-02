import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingRegressor

# 读取数据集
data = pd.read_csv("descriptors.csv")

selected_columns = ['BCUT2D_MWLOW', 'BCUT2D_MRLOW', 'Chi1n', 'Chi4n', 'Chi4v', 'Kappa3', 'SMR_VSA1', 'VSA_EState2', 'NumHDonors', 's+', 'VIP', 'Electrophilicity_index']
selected_data1 = data[selected_columns]

# 保存选定列到新的CSV文件
selected_data1.to_csv("final_descriptors.csv", index=False)

