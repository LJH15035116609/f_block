from flask import Flask, request, jsonify
import csv
import os

app = Flask(__name__)

# 定义 first_15_metals 列表
first_15_metals = [
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu'
]

@app.route('/save_to_data_csv', methods=['POST'])
def save_to_data_csv():
    data = request.json
    metal = data['metal']
    valence = data['valence']

    # 匹配 homo.csv 中的金属和价态
    with open('homo.csv', 'r') as homo_file:
        reader = csv.reader(homo_file)
        for row in reader:
            if row[0] == metal and row[1] == valence:
                # 保存匹配的第三列到 data.csv
                with open('data.csv', 'a') as data_file:
                    writer = csv.writer(data_file)
                    writer.writerow([row[2]])
                break
    return jsonify({"status": "success"})

@app.route('/save_smiles_and_generate_descriptors', methods=['POST'])
def save_smiles_and_generate_descriptors():
    data = request.json
    smiles = data['ligandSmiles']
    metal = data['metal']

    # 保存 SMILES 到 smiles.csv
    with open('smiles.csv', 'a') as smiles_file:
        writer = csv.writer(smiles_file)
        writer.writerow([smiles])

    # 生成描述符
    os.system('python all_descriptors.py')

    # 特征选择
    if metal in first_15_metals:
        os.system('python descriptors.py')

    # 绘图
    os.system('python draw.py')

    # 假设生成的图像保存为 image.png
    image_url = '/path/to/generated/image.png'  # 修改为你的图像路径
    return jsonify({"image_url": image_url})

@app.route('/predict_logk', methods=['POST'])
def predict_logk():
    data = request.json
    ionic_strength = data['ionicStrength']
    metal = data['metal']

    # 读取 data.csv，并附加 ionic_strength
    with open('data.csv', 'a') as data_file:
        writer = csv.writer(data_file)
        writer.writerow([ionic_strength])

    # 预测 logK
    if metal in first_15_metals:
        os.system('python ac_model.py --mode=first15')
    else:
        os.system('python ac_model.py --mode=last15')

    # 假设 ac_model.py 生成了一个结果文件 result.txt
    with open('result.txt', 'r') as result_file:
        prediction = result_file.read()

    return jsonify({"prediction": prediction})

if __name__ == '__main__':
    app.run(debug=True)
