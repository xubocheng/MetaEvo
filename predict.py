import pandas as pd
import torch
from transformers import AutoTokenizer
import os
import torch
from torch.utils.data import DataLoader, Dataset
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import pandas as pd
from tqdm import tqdm
import os

def seq_text(df=None):
    text_data = [x.strip() for x in df['Sequence'].tolist()]
    print(len(text_data))
    return text_data    


def spl_seq(seq=None,n=None):
    l = len(seq)
    step = int(l/n)
    print(step)
    subseq = [seq[i:i+step] for i in range(0, l, step)]
    return subseq


def mercsvfile_for_multipred(title,item=None,n=None,gpu_id = None):
    data = pd.DataFrame()
    for i in range(n):
        item=str(i)
        filename = 'intermediate_result/pos_neg_item/{}{}_{}.csv'.format(title, item, gpu_id)
        df = pd.read_csv(filename, sep=',')
        data = pd.concat([data,df])
        i+=1
    data.reset_index(drop=True, inplace=True)
    data = data.sort_values(by='log10_mic', ascending=True)
    data.to_csv('intermediate_result/round_initial/{}all_{}_{}.csv'.format(title, item, gpu_id), index=False, header=True)
    return data

class SequencesDataset(Dataset):
    """用于tokenize序列的数据集类"""
    def __init__(self, sequences, tokenizer, max_length=20):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        inputs = self.tokenizer(sequence, return_tensors="pt", max_length=self.max_length, truncation=True, padding="max_length")
        return {
            'input_ids': inputs['input_ids'],
            'attention_mask': inputs['attention_mask'],
        }

def load_and_process_model(model_path, dataloader, is_regression=False):
    """加载模型并处理数据"""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if is_regression:
        
        model = AutoModelForSequenceClassification.from_pretrained(model_path, num_labels=2 if not is_regression else 1).to(device)
        new_weights_path = './model_weights/reg_pytorch_model.bin'
        #new_weights_path = '/data/wangaw/ESM/output/broad_regression_dropout/checkpoint-4000/pytorch_model.bin'
    else:
        model = AutoModelForSequenceClassification.from_pretrained(model_path, num_labels=2 if not is_regression else 1).to(device)
        new_weights_path = './model_weights/classify_pytorch_model.bin'

    if os.path.exists(new_weights_path):
        new_state_dict = torch.load(new_weights_path)
        model.load_state_dict(new_state_dict)
        print(f"Loaded model weights from {new_weights_path}")
    else:
        print(f"Warning: Model weights file not found at {new_weights_path}. Using the original model weights from {model_path}.")
        
    results = batched_process(model, dataloader, is_regression=is_regression)
    del model  # 显式删除模型以释放显存
    torch.cuda.empty_cache()  # 清空未使用的缓存
    return results

def batched_process(model, dataloader, is_regression=False):
    """分批处理数据"""
    model.eval()
    results = []
    
    with torch.no_grad():
        for batch in tqdm(dataloader, desc="Processing batches"):
            input_ids = batch['input_ids'].squeeze(1).to(model.device)
            attention_mask = batch['attention_mask'].squeeze(1).to(model.device)
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            if is_regression:

                labels = outputs.logits.squeeze()
                # 将单个浮点数转换为列表
                if isinstance(labels, float):
                    results.extend([labels])
                else:
                    results.extend(labels.tolist())

            else:
                preds = torch.argmax(outputs.logits, dim=1).tolist()
                results.extend(preds)
    return results

def filter_sequences_with_labels(sequences, predictions, variance=None, threshold=None, is_regression=False):
    """根据预测过滤序列并保存预测值"""
    filtered_sequences, filtered_labels, filtered_variance = [], [], []
    if is_regression:
        for seq, pred in zip(sequences, predictions):
            if pred < threshold:
                filtered_sequences.append(seq)
                filtered_labels.append(pred)
    else:
        for seq, pred in zip(sequences, predictions):
            if pred == 1:
                filtered_sequences.append(seq)
                filtered_labels.append(pred)
    return filtered_sequences, filtered_labels

def save_sequences_with_labels(sequences, labels, file_path, variance=None):
    """保存序列及其标签到文件"""
    with open(file_path, 'a') as f:
        for sequence, label in zip(sequences, labels):
            f.write(f"{sequence}\t{label}\n")

def load_existing_results(file_path):
    """加载已有的结果"""
    if not os.path.exists(file_path):
        return {}
    results = {}
    with open(file_path, 'r') as f:
        for line in f:
            seq, label = line.strip().split('\t')
            results[seq] = label
    return results

def multi_predict_sequence_in_evolution(df=None, item=False, n=None, gpu_id=None):
    all_data = pd.DataFrame()
    sequences = df['Sequence'].tolist()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    MODEL_NAME_OR_PATH = os.path.join(script_dir, 'esm_model')
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME_OR_PATH, padding_side='right', use_fast=True, model_max_length=20, trust_remote_code=True)
    classification_results={}
    # 分类预测
    if sequences:
        classification_dataset = SequencesDataset(sequences, tokenizer, max_length=20)
        classification_dataloader = DataLoader(classification_dataset, batch_size=64, shuffle=False)
        new_classification_preds= load_and_process_model(MODEL_NAME_OR_PATH, classification_dataloader, is_regression=False)
        new_classification_results = {seq: pred for seq, pred in zip(sequences, new_classification_preds)}
        classification_results.update(new_classification_results)

    positive_sequences = [seq for seq, pred in classification_results.items() if pred == 1]

    regression_results = {}
    remaining_positive_sequences = [seq for seq in positive_sequences if seq not in regression_results]

    if remaining_positive_sequences:
        regression_dataset = SequencesDataset(remaining_positive_sequences, tokenizer, max_length=20)
        regression_dataloader = DataLoader(regression_dataset, batch_size=64, shuffle=False)
        new_regression_preds = load_and_process_model(MODEL_NAME_OR_PATH, regression_dataloader, is_regression=True)
        new_regression_results = {seq: pred for seq, pred in zip(remaining_positive_sequences, new_regression_preds)}
        regression_results.update(new_regression_results)

    regression_sequences = list(regression_results.keys())
    regression_labels = list(regression_results.values())
    res_df = pd.DataFrame({
        'Sequence': regression_sequences,
        'log2_mic': regression_labels  
    })
    
    all_data = pd.concat([all_data, res_df])
    negative_sequences = [seq for seq in sequences if seq not in regression_sequences]
    negative_labels = [12] * len(negative_sequences)
    neg_df = pd.DataFrame({
        'Sequence': negative_sequences,
        'log2_mic': negative_labels
    })
    all_data = pd.concat([all_data, neg_df])
    return all_data

