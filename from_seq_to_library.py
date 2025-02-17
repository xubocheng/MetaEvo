import pandas as pd
import numpy as np
from rand_gen import *
from tqdm import tqdm
import re
# 示例用法
aa_code = {'R':'A','K':'A','H':'A','A':'B','I':'B','L':'B','V':'B','F':'C','W':'C','Y':'C','N':'D','Q':'D','S':'D','T':'D'}
def list_to_dataframe(list):
    df = pd.DataFrame(list, columns=['Sequence'])
    return df

def convert_list_to_string(df):
    df['aa_code'] = df['aa_code'].apply(lambda x: ''.join(map(str, x)))
    return df

def encode_amino_acids(dataframe, aa_code_dict, sequence_column):
    encoded_sequences = []
    for sequence in dataframe[sequence_column]:
        encoded_sequence = [aa_code_dict.get(aa, 0) for aa in sequence]
        encoded_sequences.append(encoded_sequence)

    dataframe["aa_code"] = encoded_sequences
    dataframe = convert_list_to_string(dataframe)
    return dataframe

import pandas as pd
from Bio import SeqIO
import subprocess

def create_fasta_from_dataframe(dataframe, output_filename):
    with open(output_filename, "w") as fasta_file:
        for index, row in dataframe.iterrows():
            sequence = row["aa_code"]
            fasta_file.write(f">{index}\n{sequence}\n")

    
def extract_sequences_from_fasta(input_fasta, output_csv):
    sequences = []
    with open(input_fasta, 'r') as fasta_file:
        lines = fasta_file.readlines()
        for idx in range(0, len(lines), 2):
            sequence = lines[idx + 1].strip()
            sequences.append(sequence)

    df = pd.DataFrame({'aa_code': sequences})
    df.to_csv(output_csv, index=False)


from collections import Counter

def analyze_strings_in_list(input_list, target_length, threshold, prioritize_amino_acids=False):

    group_by_length = {}
    for s in input_list:
        length = len(s)
        if length not in group_by_length:
            group_by_length[length] = []
        group_by_length[length].append(s)

    new_strings_list = []
    
    
    print(f'Frequency Threshold: {threshold}')
    
    for length, group in group_by_length.items():
        if length == target_length:
            print(f"分组：字符串长度为{length}")
            
            
            lengths = [len(s) for s in group]

          
            print(f"字符串长度统计：{Counter(lengths)}")

         
            position_frequencies = [{char: round(count/len(group), 4) for char, count in Counter(s[i] for s in group).items()} for i in range(target_length)]

            for i, frequencies in enumerate(position_frequencies):
                print(f"位置 {i+1} 的字符频率统计：{frequencies}")
            
        
            new_string = ""
            for i in range(target_length):
                sorted_chars = sorted(position_frequencies[i].items(), key=lambda x: x[1], reverse=True)

                if prioritize_amino_acids and len(sorted_chars) == 2 and sorted_chars[0][1] == sorted_chars[1][1]:
                    amino_acids = ['R', 'L', 'W']
                    chosen_char = next((char for char, _ in sorted_chars if char in amino_acids), None)
                    new_string += chosen_char if chosen_char else f"({''.join([char for char, _ in sorted_chars])})"
                else:
                    cumulative_freq = 0
                    cumulative_chars = []
                    
                    for char, freq in sorted_chars:
                        cumulative_freq += freq
                        cumulative_chars.append(char)
                        if cumulative_freq >= threshold:
                            break

                    if len(cumulative_chars) == 1:
                        new_string += cumulative_chars[0]
                    else:
                        new_string += f"({''.join(cumulative_chars)})"
            
            print(f"新字符串为：{new_string}")
            new_strings_list.append(new_string)
    
    meta_sequence = '\n'.join(new_strings_list)
    
    return meta_sequence



def generate_variations(input_string):

    start = input_string.find('(')
    end = input_string.find(')', start)

    if start == -1 or end == -1:
        return [input_string]

    before = input_string[:start]
    options = input_string[start+1:end]
    after = input_string[end+1:]

    variations = []
    for char in options:
        new_string = before + char + after
        variations.extend(generate_variations(new_string))

    return variations


def generate_library(input_list):
    all_variations = []
    input_variations_dict = {}
    variations = generate_variations(input_list)
    input_variations_dict[input_list] = variations
    return variations,input_variations_dict

def from_seq_to_library_length(input_list, target_length,threshold, prioritize_amino_acids=False):#此处放入的是cd_hit_delete_dissimilarity_seqs(df)的返回值
    new_strings_list = analyze_strings_in_list(input_list, target_length,threshold=threshold, prioritize_amino_acids=prioritize_amino_acids)#获取meta peptides
    
    all_variations,new_strings_list_variations_dict = generate_library(new_strings_list) #根据meta peptides获取全文库
    
    return new_strings_list, all_variations,new_strings_list_variations_dict

def from_seq_to_library(df,l,threshold, prioritize_amino_acids=False):

    result = df['Sequence'].tolist()
    library = []
    summary_meta_seq = []
    metaseq_library_dict = {}
    meta_seq,sublibrary,seq_sublibrary_dict = from_seq_to_library_length(result, l,threshold=threshold, prioritize_amino_acids=prioritize_amino_acids)
    
    library.extend(sublibrary)
    metaseq_library_dict.update(seq_sublibrary_dict)
    
    return meta_seq,library,metaseq_library_dict

def generate_all_variations(input_string):
    # 定义字符映射
    character_mapping = {
        'A': ['K', 'R', 'H'],
        'B': ['A', 'I', 'L', 'V'],
        'C': ['F', 'W', 'Y'],
        'D': ['N', 'Q', 'S', 'T']
    }
    all_variations = ['']

    for char in input_string:
        if char in character_mapping:
            new_variations = []
            for mapped_char in character_mapping[char]:
                for variation in all_variations:
                    new_variations.append(variation + mapped_char)

            all_variations = new_variations

    return all_variations

def generate_variations_for_list(input_list):
    all_combined_variations = []
    for item in input_list:
        all_variations = generate_all_variations(item)
        all_combined_variations.extend(all_variations)
    return all_combined_variations


def peptide_to_plibrary(df,subdf, length, threshold, prioritize_amino_acids=False):
    summary_meta_seq, code_library, meta_seq_library_dict = from_seq_to_library(subdf, length, threshold=threshold, prioritize_amino_acids=prioritize_amino_acids) # 编码的文库
    
    sub_df = df[df['Sequence'].isin(code_library)]
    mean_mic = sub_df['log2_mic'].mean()

    # 打印相关信息
    count = len(code_library)
    print(f'summary_meta_to_plibrary: {count}')
    print('mean_mic:', mean_mic)

    return summary_meta_seq, code_library, meta_seq_library_dict,sub_df, mean_mic


def merge_list(list_of_lists):
    merged_list = []
    for sublist in list_of_lists:
        merged_list.extend(sublist)
    return merged_list

def plibrary_select_subseq(meta_to_plibrary=None,H=None,limit=None):
    data = pd.DataFrame(meta_to_plibrary, columns=['Sequence'])
    df2 = pd.DataFrame(columns=data.columns)
    count = 0 
    for i in tqdm(range(6,14)):#没有必要写成14吧
        for j in tqdm(range(2,6)):
            subdf = select_seq(data=data,seq_length=i,RKH=j,limit=limit,H=H)
            df2 = df2.append(subdf, ignore_index=True)
            print(len(subdf))
            count+=len(subdf)
            print(count)
        if count >=100000:
            break
    sequences = df2.Sequence.to_list()
    return sequences

def new_plibrary_select_subseq(meta_seq_variations_dict=None,H=None,limit=None):
    filter_meta_seq_variations = {}
    for meta_seq in tqdm(meta_seq_variations_dict.keys()):
        seq_length = len(re.sub(r'\([^()]*\)', 'X', meta_seq))
        subdf_list  = []
        for j in range(2,6):
            data = pd.DataFrame(meta_seq_variations_dict[meta_seq], columns=['Sequence'])
            subdf = select_seq(data=data,seq_length=seq_length,RKH=j,limit=limit,H=H)
            if not subdf.empty:
                subdf_list.extend(subdf['Sequence'].tolist())
                
        print(f'subdf_list length',len(subdf_list))
        filter_meta_seq_variations[meta_seq] = subdf_list

    return filter_meta_seq_variations
