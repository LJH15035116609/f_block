import joblib
import pandas as pd

# 读取数据
df = pd.read_csv('data.csv')

# 加载模型
model = joblib.load('model_ac_logkl.pkl') 

predictions= model.predict(df)

print(predictions)
