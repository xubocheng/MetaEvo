import pandas as pd
import numpy as np
import re
from rand_gen import *
from from_seq_to_library import *
import warnings
warnings.filterwarnings("ignore")

import sys
import os
from predict import *

def first_pipeline(num,net,l):
    df = list_to_dataframe(merge_meta(num,net,l))
    
    sequence = df['Sequence'].to_list()
    
    num_random_sequence = len(sequence)

    print('共生成meta序列:',num_random_sequence)

    return sequence

def process_sequences(file_path):
    with open(file_path, 'r') as file:
        sequences = file.readlines()
    
    for i, sequence in enumerate(sequences):
        modified_sequence = generate_variations(sequence.strip())
        filename = f'./data/library/random_metaseq_library_{i}.txt'
        with open(filename, 'w') as file:
            file.write('Sequence\n')  
            modified_sequence_str = '\n'.join(modified_sequence)
            file.write(modified_sequence_str + '\n')

    return print(f'共处理序列数量：{len(sequences)}')

def meta_predict_evolution(meta_seq,gpu_id):

    res_df,broad_AMPs = predict_meta_mic(meta_seq,gpu_id)
    result_df = meta_evolution(meta_seq,res_df,broad_AMPs,gpu_id)
    
    return result_df

def predict_meta_mic(meta_seq,gpu_id):
    
    library_sequence = generate_variations(meta_seq.strip())
    num_library_sequence = len(library_sequence)
    
    
    df = pd.DataFrame(library_sequence, columns=["Sequence"])
    res_df = multi_predict_sequence_in_evolution(df=df,item=False,n=1,gpu_id=gpu_id)  
    
    res_df = res_df.sort_values(by='log2_mic')

    mic_average = res_df['log2_mic'].mean()
    broad_AMPs = res_df[res_df['log2_mic']<=7]
    num_broad_AMPs = broad_AMPs.shape[0]
    
    print('meta序列：', meta_seq)
    print('meta序列文库数量：', num_library_sequence)
    print('meta序列log10_mic：', mic_average)    
    print('筛选广谱性AMPs个数：', num_broad_AMPs) 
    return res_df,broad_AMPs

def meta_evolution(meta_seq,res_df,broad_AMPs,gpu_id):
    num_broad_AMPs = broad_AMPs.shape[0]
    #num = len(res_df)
   
    mic_average = res_df['log2_mic'].mean()
    print(meta_seq,'开始进化')
    print('初始log10_mic：',mic_average)

    data = {
        "随机肽参数": {
        "meta序列": meta_seq,
        "初始log2_mic": mic_average,
    },
        "预测与筛选结果": [{
            "进化轮次":'0',
            "抗菌肽鸡尾酒序列":meta_seq,
            "文库数量":res_df.shape[0],
            "log2_MIC":mic_average}            
        ]
    }
    
    sub_df = broad_AMPs
    
    for i in range(1,11):
        print(f'第{i}轮')
        sub_df = sub_df.head(int(num_broad_AMPs))
        
        summary_meta_seq, code_library, meta_seq_library_dict, sub_df, mean_mic = peptide_to_plibrary(
            res_df, sub_df, length=13, threshold=0.7, prioritize_amino_acids=True
        )
        
        
        sub_num_broad_AMPs = sub_df[sub_df['log2_mic'] <= 7].shape[0]
        num_broad_AMPs = min(num_broad_AMPs, sub_num_broad_AMPs)
        if num_broad_AMPs >= 50000:
            num_broad_AMPs *= 0.05
        elif num_broad_AMPs >= 20000:
            num_broad_AMPs *= 0.1
        elif num_broad_AMPs >= 10000:
            num_broad_AMPs *= 0.2
        elif num_broad_AMPs >= 5000:
            num_broad_AMPs *= 0.3
        else:
            num_broad_AMPs *= 0.5

        
        record_txt = {
                "进化轮次":i,
                "抗菌肽鸡尾酒序列":summary_meta_seq,
                "文库数量":len(code_library),
                "log2_MIC":mean_mic
        }
        data["预测与筛选结果"].append(record_txt)
        print(f"进化轮次: {i}")
        print(f"Meta序列: {summary_meta_seq}")
        print(f"文库数量: {len(code_library)}")
        print(f"log2_MIC: {mean_mic}")
        if len(code_library) <=4:
            break       
            

    
    current_dir = os.path.dirname(os.path.abspath(__file__))

    result_folder_name = "esm_result"
    result_folder_path = os.path.join(current_dir, result_folder_name)

    if not os.path.exists(result_folder_path):
        os.makedirs(result_folder_path)

    file_name = f"{meta_seq}_evolution_debug.txt"
    file_path = os.path.join(result_folder_path, file_name)
    save_and_print_data(data, file_path)    
    

    final_df = out_to_df(data)
    print('完成',meta_seq,'进化')
    return final_df

def save_and_print_data(data, file_path):
    with open(file_path, 'w', encoding='utf-8') as file:
        for round_data in data["预测与筛选结果"]:
            content = (
                f"  进化轮次: {round_data['进化轮次']}\n"
                f"  Meta序列: {round_data['抗菌肽鸡尾酒序列']}\n"
                f"  文库数量: {round_data['文库数量']}\n"
                f"  log2_MIC: {round_data['log2_MIC']}\n\n"
            )
            # 写入文件
            file.write(content)
            # 打印到控制台
            print(content, end='')        
            

def extract_first_column(file_path):
        df = pd.read_csv(file_path)
        first_column_list = df['initial_meta'].tolist()
        return first_column_list

def extract_existing_sequences(folder_path):
    existing_sequences = []
    for filename in os.listdir(folder_path):
        if filename.endswith('_evolution_debug.txt'):
            sequence = filename.split('_evolution_debug.txt')[0]
            existing_sequences.append(sequence)
    return existing_sequences



def extract_and_sort(folder_path):
    # 空的数据框，用于存储文件名和数字
    df = pd.DataFrame(columns=['Meta_sequence', 'log10_MIC'])

    # 遍历指定文件夹下的所有文件
    for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)

            # 读取文件内容
            with open(file_path, 'r') as file:
                lines = file.readlines()
                if lines:
                    # 提取最后一行中的数字（浮点数）
                    last_line = lines[-2]
                    numbers = re.findall(r'\d+\.\d+|\d+', last_line)
                    if numbers:
                        number = float(numbers[-1])  
                        df = df.append({'Meta_sequence': filename.replace('_evolution.txt', ''), 'log10_MIC': number}, ignore_index=True)

    # 对数字列进行排序
    df = df.sort_values(by='log10_MIC')
    return df

#汇总结果
def summary_txt_result():
    path = os.getcwd()
    folder_path = path+'/final_result/'
    result_df = extract_and_sort(folder_path)
    result_df.to_csv(folder_path+'summary_result_{}.csv'.format(result_df.shape[0]),index=False,header=True)
    print('汇总结果已保存')
    return result_df

def get_existing_meta_sequences(folder_path):
    existing_meta_sequences = []
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            meta_sequence = file.replace('.csv', '')
            existing_meta_sequences.append(meta_sequence)
    return existing_meta_sequences

def out_to_df(data):
    rounds = data["预测与筛选结果"]
    initial_data = rounds[0]
    initial_meta = initial_data['抗菌肽鸡尾酒序列']
    initial_library_num = initial_data['文库数量']

    last_round = rounds[-1]
    if last_round['文库数量'] < 4:
        final_round_data = rounds[-2] if len(rounds) > 1 else rounds[-1]  # 兼容只有一轮的情况
    else:
        final_round_data = last_round
        
    final_meta = final_round_data['抗菌肽鸡尾酒序列']
    final_library_num = final_round_data['文库数量']
    final_mic = final_round_data['log2_MIC']
    
    # 创建包含单行数据的DataFrame
    df = pd.DataFrame({
        'initial_meta': [initial_meta],
        'initial_library_num': [initial_library_num],
        'final_meta': [final_meta],
        'final_library_num': [final_library_num],
        'log2_MIC': [final_mic]
    })
    
    return df
