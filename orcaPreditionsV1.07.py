#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Script to predict 3D organization from DNA sequence using ORCA.
"""

import sys 
import pandas as pd
import os
import numpy as np
import gzip
import time
import argparse
import logging
from Bio import SeqIO
# Import required libraries.
# If they are not available, propagate ImportError with a helpful message for Genotoul users.
try:
    from selene_sdk.sequences import Genome
    import orca_predict
except ImportError as e:
    sys.stderr.write(
        "Error: Missing required packages ('selene_sdk' or 'orca_predict').\n"
        "To run this script on Genotoul, make sure to load the required modules:\n"
        "    module load devel/Miniconda/Miniconda3\n"
        "    module load bioinfo/Orca/88df3b5\n\n"
        "If you still encounter this error, please run the script using the environment's Python executable directly:\n"
        "    /usr/local/bioinfo/src/Miniconda/Miniconda3/envs/orca-88df3b5_env/bin/python orcaPreditionsV1.07.py [args]\n\n"
    )
    raise e

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

def parse_args():
    parser = argparse.ArgumentParser(description="Predict 3D organization from DNA sequence using ORCA.", epilog="""
# Input tsv file format:
------------------------
Tabulated file with header, containing as many lines as there are predictions to be made.
Mandatory columns are:
    ID: unique ID which must match fasta sequence file name: E.g, if ID = 'ID1', sequence fasta file must be 'ID1.fa' (or 'ID1.fa.gz' and 'ID1.fa.gzip' for 1M model).
    model: 32000000 or 1000000 ORCA model (i.e., 32M or 1M respectively).
    scale: list of prediction scales to keep. 
        For all 32M model scales, use: '32000000,16000000,8000000,4000000,2000000,1000000'.
        For 1M model, scale must be 1000000.
    wpos: The coordinate of the center position of the sequence, 
        which is start position + 16000000 or + 500000 for 32M and 1M model respectively.
    mpos: The coordinate to zoom into for multiscale prediction of 32M model. 
        mpos = 16000000 (for 32M model) to zoom at the center. Not used for 1M model but mandatory.
    model_HFF: 0 or 1 to save predictions from HFF model.
    model_ESC: 0 or 1 to save predictions from ESC model.

# Versions:
-----------
    1.00:   predictions with 1M and 32M Orca models on CPU mode.
    1.02:   add GPU mode.
    1.03:   Check columns values for ['model_HFF', 'model_ESC', 'model']
            Bug fix: when wpos != mpos
    1.04:   possibility to read compressed fasta files with extentions: .fa.gz and .fa.gzip (for 1M model)
    1.05:   possibility to compress (.gz) predicted matrices (for 1M model)
    1.06:   read and write compressed fasta and matrices respectively for 32M model
    1.07:   add option to write log files in real time, add argparse, optimize I/O with pandas
    Next:
        -Check the name of the Fasta sequence “>ID,” because if it doesn't exist, there's a bug
        -32M: add the option of keeping only 'predictions' (i.e Obs/Exp) and/or 'normmat' (i.e Obs) matrices,
        -add 256Mb model (?),
        -add default parameters for missing input.tsv columns (scales, mpos, model_HFF, model_ESC...).

# Examples:
-----------
    python orcaPredictionsV1.07.py -i input.tsv -f seq_folder -o output_folder
    
    ## CPU mode on genotoul:
    
    srun --mem 32G --pty bash
    module load devel/Miniconda/Miniconda3
    module load bioinfo/Orca/88df3b5
    /work/user/nmary/script/orcaPredictionsV1.07.py \\
        --input /work/user/nmary/script/Datas_example/OrcaPredictions/control.tsv \\
        --fasta_dir /work/user/nmary/script/Datas_example/OrcaPredictions/Sequences \\
        --output ~/work/tmp
        
    ## GPU mode on genotoul:
    
    echo '#!/bin/bash
    #SBATCH -p gpuq
    #SBATCH -J orca 
    #SBATCH -o orca.out 
    #SBATCH -e orca.err 
    #SBATCH -t 00-00:50:00 # "days-hours:minutes:seconds"
    #SBATCH --gres=gpu:A100_2g.20gb:1
    #SBATCH --mem=32G

    #Load modules
    #Need Miniconda
    module load devel/Miniconda/Miniconda3
    module load bioinfo/Orca/88df3b5

    /work/user/nmary/script/orcaPredictionsV1.07.py \\
        --input /work/user/nmary/script/Datas_example/OrcaPredictions/control.tsv \\
        --fasta_dir /work/user/nmary/script/Datas_example/OrcaPredictions/Sequences \\
        --output ~/work/tmpGPU \\
        --cuda' > orca.sh

    chmod u+x orca.sh
    sbatch orca.sh
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--input", "-i", type=str, required=True,
                        help="Tab-delimited file with prediction metadatas (input TSV file).")
    parser.add_argument("--fasta_dir", "-f", type=str, required=True,
                        help="Directory containing fasta sequence files.")
    parser.add_argument("--output_dir", "-o", type=str, required=True,
                        help="Directory to save output files (will be created if required).")
    parser.add_argument("--cuda", action="store_true",
                        help="Enable CUDA GPU.")
    parser.add_argument("--compress", action="store_true",
                        help="Compress predicted matrices (.gz).")
    
    return parser.parse_args()


def validate_inputs(args):
    input_tsv = os.path.expanduser(args.input)
    seq_folder = os.path.expanduser(args.fasta_dir)
    output_folder = os.path.expanduser(args.output_dir)
    
    if not os.path.isfile(input_tsv):
        logging.error(f"Input TSV file '{input_tsv}' not found.")
        exit(1)
        
    if not os.path.isdir(seq_folder):
        logging.error(f"Sequence folder '{seq_folder}' not found.")
        exit(1)
    
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    # Checking input_tsv
    required_columns = {'ID', 'model', 'scale', 'wpos', 'mpos', 'model_HFF', 'model_ESC'}
    try:
        df = pd.read_csv(input_tsv, sep='\t')
    except Exception as e:
        logging.error(f"Failed to read input TSV: {e}")
        exit(1)

    # Checking mandatory columns in input_tsv    
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        logging.error(f"Missing required columns in input TSV: {', '.join(missing_columns)}")
        exit(1)
    
    # Check whether columns ('model_HFF', 'model_ESC') contain only 1 or 0 :    
    for column in ['model_HFF', 'model_ESC']:
        invalid_rows = df[~df[column].isin([0, 1])]
        if not invalid_rows.empty:
            logging.error(f"Column '{column}' contains values other than 0 or 1. IDs: {', '.join(map(str, invalid_rows['ID'].tolist()))}")
            exit(1)
    
    # Check whether column ('model') contain only 32000000 or 1000000 :
    for column in ['model']:
        invalid_rows = df[~df[column].isin([32000000, 1000000])]
        if not invalid_rows.empty:
            logging.error(f"Column '{column}' contains values other than 32000000 or 1000000. IDs: {', '.join(map(str, invalid_rows['ID'].tolist()))}")
            exit(1)
    
    df = df.sort_values(by='model', ascending=False).reset_index(drop=True)
    return df, seq_folder, output_folder

#print + reset timer
start_time = time.time()
def log_and_reset_timer(message):
    """
    add log with elapsed time between each function call
    message = beginning of sentence
    """
    global start_time
    elapsed_time = time.time() - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.info(f"{message} {int(hours)}h {int(minutes)}m {int(seconds)}s")
    # reset timer
    start_time = time.time()


model_dict = {32000000: '32M', 1000000: '1M'}
def load_orca_model(model, use_cuda):
    """Load ORCA model and return OrcaModel (i.e. '32M' or '1M')."""
    OrcaModel = model_dict[model]
    logging.info(f"Loading {OrcaModel} model...")
    orca_predict.load_resources(models=OrcaModel, use_cuda=use_cuda)
    log_and_reset_timer("Model loaded successfully in")
    return OrcaModel


fasta_extensions = ['.fa', '.fa.gz', '.fa.gzip']
def get_sequence(ID, seq_dir):
    """Reads and returns the DNA sequence."""
    fasta_prefix = os.path.join(seq_dir, str(ID))
    for ext in fasta_extensions:
        fasta_path = fasta_prefix + ext
        if os.path.isfile(fasta_path):
            if fasta_path.endswith(('.fa.gz', '.fa.gzip')):
                with gzip.open(fasta_path, "rt") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                    if not records:
                        logging.error(f"No sequences found in '{fasta_path}'.")
                        exit(1)
                    return str(records[0].seq)
            else:
                records = list(SeqIO.parse(fasta_path, "fasta"))
                if not records:
                    logging.error(f"No sequences found in '{fasta_path}'.")
                    exit(1)
                return str(records[0].seq)
    
    logging.error(f"No sequences found for '{fasta_prefix}' in directory.")
    exit(1)


scale_dict = {'32000000': 0, '16000000': 1, '8000000': 2, '4000000': 3, '2000000': 4, '1000000': 5}

# export matrix with panda
def export_matrix_fast(matrix, filepath, compressed):
    df_mat = pd.DataFrame(matrix)
    if compressed:
        df_mat.to_csv(filepath + '.gz', sep='\t', header=False, index=False, compression='gzip')
    else:
        df_mat.to_csv(filepath, sep='\t', header=False, index=False)


def export32M_matrix(row, OrcaModel, orcaCellModel, outputs, compressed_matrix, output_dir):
    """
    From OrcaModel 32M predictions outputs:
    export normat and prediction matrices for scale matrices to keep (as .tsv)
    return scale_indexes (output indexes to keep) use in export32M_parameters function
    """
    scales = str(row['scale']).split(',')
    # check scales values
    try:
        scale_indexes = [scale_dict[scale] for scale in scales]
    except KeyError as e:
        logging.error(f"{row['ID']} Error: Invalid scale '{e.args[0]}' in the input TSV file.")
        exit(1)
        
    for index in scale_indexes:
        #Exports prediction and normat matrices
        matrix_pred = outputs['predictions'][0][index]
        matrix_norm = outputs['normmats'][0][index]
        suffix = list(scale_dict.keys())[index] + '_' + OrcaModel + '_' + orcaCellModel + '.tsv'

        # Save matrix as .tsv files
        pred_file_path = os.path.join(output_dir, f"{row['ID']}_predictions_{suffix}")
        norm_file_path = os.path.join(output_dir, f"{row['ID']}_normmats_{suffix}")
        
        export_matrix_fast(matrix_pred, pred_file_path, compressed_matrix)
        export_matrix_fast(matrix_norm, norm_file_path, compressed_matrix)
        
    logging.info(f"{OrcaModel}_{orcaCellModel} predictions exported successfully for {row['ID']}.")
    return scale_indexes


def export32M_parameters(row, OrcaModel, outputs, scale_indexes, output_dir):
    """
    From OrcaModel 32M predictions outputs:
    export matrix parameters (start, stop...) for scale matrices to keep (as .bed).
    """
    param_dict = {'ID': [], 'start': [], 'end': [], 'bin_width': [], 'mpos': [], 'wpos': [], 'model': []}
    for index in scale_indexes:
        param_dict['ID'].append(str(row['ID']) + '_predictions_' + list(scale_dict.keys())[index])
        param_dict['start'].append(outputs['start_coords'][index])
        param_dict['end'].append(outputs['end_coords'][index])
        param_dict['bin_width'].append(str((outputs['end_coords'][index] - outputs['start_coords'][index]) / 250))
        param_dict['mpos'].append(row['mpos'])
        param_dict['wpos'].append(row['wpos'])
        param_dict['model'].append(row['model'])
    
    path_bed = os.path.join(output_dir, f"{row['ID']}.bed")
    pd.DataFrame.from_dict(param_dict).to_csv(path_bed, header=True, sep='\t', index=False)
    logging.info(f"{OrcaModel} parameters exported successfully for {row['ID']}.")


def export1M_matrix_parameters(row, OrcaModel, orcaCellModel, outputs, compressed_matrix, output_dir):
    """
    From OrcaModel 1M predictions outputs:
    export matrix and parameters (start, stop...) as .bed
    """
    if not str(row['scale']) == '1000000':
        logging.error(f"{row['ID']} Error: Invalid scale '{row['scale']}' in the input TSV.")
        exit(1)
    
    suffix = str(row['ID']) + '_predictions_' + str(row['scale']) + '_' + OrcaModel + '_' + orcaCellModel + '.tsv'
    pred_file_path = os.path.join(output_dir, suffix)

    export_matrix_fast(outputs, pred_file_path, compressed_matrix)
    
    param_dict = {
      'ID': [str(row['ID']) + '_predictions_1000000'],
      'start': [row['wpos'] - 500000],
      'end': [row['wpos'] + 500000],
      'bin_width': [4000.0],
      'mpos': [row['mpos']],
      'wpos': [row['wpos']],
      'model': [row['model']]
    }
    path_bed = os.path.join(output_dir, f"{row['ID']}.bed")
    pd.DataFrame.from_dict(param_dict).to_csv(path_bed, header=True, sep='\t', index=False)
    logging.info(f"{OrcaModel}_{orcaCellModel} parameters exported successfully for {row['ID']}.")


def process_predictions(input_tsv_df, seq_dir, output_dir, use_cuda, compressed_matrix):
    """Processes predictions for each row of the input TSV."""
    OrcaModel = load_orca_model(input_tsv_df['model'][0], use_cuda)
    
    for i, row in input_tsv_df.iterrows():
        # Change ordered OrcaModel if needed
        if OrcaModel != model_dict[row['model']]:
            #end timer for the first OrcaModel prediction batch
            log_and_reset_timer(f"{OrcaModel} ORCA predictions successfully completed in")
            OrcaModel = load_orca_model(row['model'], use_cuda)
        
        logging.info(f"Processing {row['ID']}...")
        
        # Validate sequence length
        sequence = get_sequence(row['ID'], seq_dir)
        if len(sequence) != row['model']:
            logging.error(f"Sequence length ({len(sequence)}) does not match model size ({row['model']}) for {row['ID']}.")
            exit(1)
        
        encoded_sequence = Genome.sequence_to_encoding(sequence)[None, :, :]
        
        ################## 32M #########################
        if OrcaModel == '32M':
            # ORCA prediction hff
            if row['model_HFF'] == 1:
                outputs = orca_predict.genomepredict(
                    mchr=row['ID'], mpos=row['mpos'], wpos=row['wpos'],
                    models=['hff'], sequence=encoded_sequence, use_cuda=use_cuda
                )
                # Export matrices hff
                scale_idx = export32M_matrix(row, OrcaModel, 'hff', outputs, compressed_matrix, output_dir)
            
            # ORCA prediction h1esc
            if row['model_ESC'] == 1:
                outputs = orca_predict.genomepredict(
                    mchr=row['ID'], mpos=row['mpos'], wpos=row['wpos'],
                    models=['h1esc'], sequence=encoded_sequence, use_cuda=use_cuda
                )
                # Export matrices h1esc
                scale_idx = export32M_matrix(row, OrcaModel, 'esc', outputs, compressed_matrix, output_dir)
            
            #export parameters (same parameters for esc or hff)
            export32M_parameters(row, OrcaModel, outputs, scale_idx, output_dir)
        
        elif OrcaModel == '1M':
            tensor_seq = orca_predict.torch.FloatTensor(encoded_sequence).transpose(1, 2)
            if use_cuda:
                tensor_seq = tensor_seq.cuda()
                
            if row['model_HFF'] == 1:
                outputs = orca_predict.hff_1m(tensor_seq).squeeze().detach().cpu().numpy()
                export1M_matrix_parameters(row, OrcaModel, 'hff', outputs, compressed_matrix, output_dir)
                
            if row['model_ESC'] == 1:
                outputs = orca_predict.h1esc_1m(tensor_seq).squeeze().detach().cpu().numpy()
                export1M_matrix_parameters(row, OrcaModel, 'esc', outputs, compressed_matrix, output_dir)
        
    log_and_reset_timer(f"{OrcaModel} ORCA predictions successfully completed in")


if __name__ == "__main__":
    args = parse_args()
    input_tsv_df, seq_dir, output_dir = validate_inputs(args)
    process_predictions(input_tsv_df, seq_dir, output_dir, args.cuda, args.compress)
