import argparse
from metaseq_peptide import meta_predict_evolution, summary_txt_result,first_pipeline
import os
import re
import pandas as pd
def main(meta_sequence):
    gpu_id = os.environ.get('CUDA_VISIBLE_DEVICES', 'default_gpu_id')
    result_df = meta_predict_evolution(meta_sequence, gpu_id=gpu_id)
    # The format is as follows:
    # #df = pd.DataFrame({
    #     'initial_meta': [initial_meta],
    #     'initial_library_num': [initial_library_num],
    #     'final_meta': [final_meta],
    #     'final_library_num': [final_library_num],
    #     'log2_MIC': [final_mic]
    # })
    print(f"Prediction of meta sequence '{meta_sequence}' completed.")
    return result_df
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run metaseq peptide predictions.")
    parser.add_argument("--meta_sequence",
                        type=str,
                        default=None,
                        help="The meta sequence to be predicted.")
    parser.add_argument("--cationic_amino_acid_number",
                        type=int,
                        default=None,
                        help="The number of cationic amino acids for generating random meta sequences.")

    args = parser.parse_args()
    if args.meta_sequence:
        # Case 1: User provides meta_sequence
        main(args.meta_sequence)
        
    elif args.cationic_amino_acid_number:
        # Case 2: User does not provide meta_sequence, but provides cationic_amino_acid_number
        random_meta_sequences = first_pipeline(num=1, net=args.cationic_amino_acid_number, l=13)
        
        all_results = []
        for meta_seq in random_meta_sequences:
            result_df = main(meta_seq)
            all_results.append(result_df)
        
        # You can choose to merge all results into one DataFrame
        if all_results:
            combined_results_df = pd.concat(all_results, ignore_index=True)
            print("The prediction results of all randomly generated meta sequences have been merged.")
            print(combined_results_df)

    else:
        print("Error: You must provide --meta_sequence or --cationic_amino_acid_number.")
        parser.print_help()
# Example: (QTS)(KR)(KR)(ILV)(ILV)(KR)(ILV)(KR)(ILV)(ILV)(KR)(ILV)(ILV)