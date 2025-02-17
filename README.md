
# MetaEvo Project

AI-driven approaches have significantly accelerated the discovery of antimicrobial peptides, which are widely regarded as next-generation antibiotics. Higher-order antimicrobial peptide cocktails (AMPCs) that mimic host-derived combinatorial mechanisms to address monotherapy resistance risks and limitations offer a promising strategy. However, de novo design of AMPCs is an unmet challenge. Here, we present MetaEvo, an intelligent design platform for AMPCs that integrates a sequence generator, plug-and-play AI modules, and evolutionary algorithms, enabling the sustainable discovery and auto-directed evolution of AMPCs. MetaEvo achieves outstanding performance on hit rates of highly active AMPCs, module self-upgradability, and data augmentation. Using MetaEvo, we generated 11,993 templates from a 14,725,250,460-sequence space, identified 4,111 candidate AMPCs, synthesized and empirically tested 225 mixtures, and ultimately identified 105 AMPCs. The lead AMPCs exhibit potent antibiotic activity against clinically relevant drug-resistant pathogens and anti-infective efficacy in mouse models of lethal A. baumannii peritonitis and cecal ligation and puncture. MetaEvo enables the sustainable discovery of cost-effective AMPCs, provides insights into high-throughput AMP screening and the exponential proliferation of AMP repertoire, and presages a coming era in AMPs translation.

## Overview

- **Random Meta-Peptide Sequence Generation**: Generates a specified number of random Meta-peptide sequences.
- **Plug-and-Play machine learning Prediction**: Predicts and filters the sequences by finetuned protein language models.
- **Rule-based Automated Sequence Evolution**: generate AMP cocktails with controllable library sizes and high activity from vast, disordered sequence spaces


## Environment Setup

1. **Install Conda**

   Ensure that you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) installed.

2. **Create a Conda Environment**

   Open a terminal and run the following command to create a new Conda environment:

   ```bash
   conda create -n my_project_env python=3.8
   ```

   This will create a Conda environment named `my_project_env` and install Python 3.8.

3. **Activate the Conda Environment**

   Use the following command in the terminal to activate your new environment:

   ```bash
   conda activate my_project_env
   ```
4. **Install Required Libraries via pip**

   While in the activated environment, run the following command to install all required libraries:

   ```bash
   pip install pandas numpy tqdm biopython transformers==4.33.0 accelerate==0.22.0 torch==2.0.1
   ```

5. **Download Model Files**

   - Download the [ESM-2](https://huggingface.co/facebook/esm2_t33_650M_UR50D) model files into an `esm_model` folder located in the project directory.
   - Download the [classification model](https://zenodo.org/records/14880503) and [regression model](https://zenodo.org/records/14880503) weights from Zenodo and save them in a `model_weights` folder in the project directory.

## Usage

MetaEvo is executed through `meta_user.py`, providing two different input modes:

1. **Random Generation Mode**

   ```bash
   CUDA_VISIBLE_DEVICES=0 python meta_user.py --cationic_amino_acid_number=5
   ```

   In this mode, the system will randomly generate 8 Meta sequences with the cationic amino acids number set to 5, and perform prediction and evolution on them.

2. **User-Provided Sequence Mode**

   ```bash
   CUDA_VISIBLE_DEVICES=0 python meta_user.py --meta_sequence '(QTS)(KR)(KR)(ILV)(ILV)(KR)(ILV)(KR)(ILV)(ILV)(KR)(ILV)(ILV)'
   ```

   In this mode, users can directly provide a Meta sequence through the command line, and the system will use this sequence as input for prediction and evolution.

## Code Modules

### 1. `meta_user.py`

- Handles command line inputs and calls the prediction and evolution functions.
- Uses one of the two input methods to determine whether to call `first_pipeline` for generating random sequences or directly call `meta_predict_evolution`.

### 2. Core Evolution Process

#### a. `meta_predict_evolution`

- **Input**: User-specified Meta sequence or generated random sequences.
- **Function**: Integrates `predict_meta_mic` and `meta_evolution` modules.
- **Output**: Predicted results and evolved sequences.

#### b. `predict_meta_mic`

- **Function**: Conducts classification and regression predictions on the Meta library.

#### c. `meta_evolution`

- **Function**: Conducts rule-based automatic evolution to optimize peptide sequence performance.
- Performs multiple evolution rounds to continuously filter and enhance antimicrobial peptide sequences.





