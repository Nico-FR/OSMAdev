#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Script to predict 3D organization from DNA sequence using ORCA.
"""

import sys 
import pandas as pd
import os
import numpy as np
from selene_sdk.sequences import Genome
import orca_predict
#from orca_predict import *
from Bio import SeqIO
import gzip
import time

def print_help_and_exit():
    """Displays help information and exits."""
    help_text = """
    Usage: orcaPredictions.py [input tsv file] [input sequence folder] [output prediction folder] [bool Use CUDA GPU] [bool gzip output matrices]

    # Parameters:
    -----------
    1- input tsv file: Tab-delimited file with prediction parameters (see details below).
    2- input sequence folder: Directory containing fasta sequence files.
    3- output prediction folder: Directory to save output files (will be created if required).
    4- Use CUDA GPU: Boolean (True/False) to enable CUDA GPU.
    5- Compress predicted matrices: Boolean (True/False) to compress matrices.
    
    # Input tsv file:
    ----------
    Tabulated file with header, containing as many lines as there are predictions to be made.
        Mandatory columns are:
            ID: unique ID which must match fasta sequence file: E.g, if ID = 'ID1', sequence fasta file must be = 'ID1.fa' and sequence name = ">ID1".
            model: 32000000 or 1000000 ORCA model (i.e., 32M or 1M respectively).
            scale: list of prediction scales to keep. 
                For all 32M model scales, use: '32000000,16000000,8000000,4000000,2000000,1000000'.
                For 1M model, scale must be 1000000.
            wpos: The coordinate of the center position of the sequence, 
                which is start position + 16000000 or + 5000000 for 32M and 1M model respectively.
            mpos: The coordinate to zoom into for multiscale prediction of 32M model. 
                mpos = 16000000 (for 32M model) to zoom at the center. Not used for 1M model but mandatory column.
            model_HFF: 0 or 1 to save predictions from HFF model.
            model_ESC: 0 or 1 to save predictions from ESC model.
    
    # Versions:
    --------
    1.00:   predictions with 1M and 32M Orca models on CPU mode.
    1.02:   add GPU mode.
    1.03:   Check columns values for ['model_HFF', 'model_ESC', 'model']
            Bug fix: when wpos != mpos
    1.04:   possibility to read compressed fasta files with extentions: .fa.gz and .fa.gzip (for 1M model)
    1.05:   possibility to compress (.gz) predicted matrices (for 1M model)
    1.06:   read and write compressed fasta and matrices respectively for 32M model
    Next:
        -véridier le nom de la sequance fasta '>ID' car si elle n'existe pas il y a un bug
        -32M: add the option of keeping only 'predictions' (i.e Obs/Exp) and/or 'normmat' (i.e Obs) matrices,
        -add 256Mb model (?),
        -add default parameters for missing input.tsv columns (scales, mpos, model_HFF, model_ESC...).
    
    # Examples:
    --------
    orcaPredictions.py input.tsv seq_folder output_folder True False
    
    ## CPU mode on genotoul:
    
    srun --mem 32G --pty bash
    module load devel/Miniconda/Miniconda3
    module load bioinfo/Orca/88df3b5
    /work/user/nmary/script/orcaPreditionsV1.02.py \\
        /work/user/nmary/script/Datas_example/OrcaPredictions/control.tsv \\
        /work/user/nmary/script/Datas_example/OrcaPredictions/Sequences \\
        ~/work/tmp False False
        
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

    /work/user/nmary/script/orcaPreditionsV1.02.py /work/user/nmary/script/Datas_example/OrcaPredictions/control.tsv /work/user/nmary/script/Datas_example/OrcaPredictions/Sequences ~/work/tmpGPU True False' > orca.sh

    chmod u+x orca.sh
    sbatch orca.sh
    """
    print(help_text)
    exit(1)

def validate_inputs():
    """
    Validates command-line arguments and read input files.
    Return all arguments
    """
    if len(sys.argv) != 6:
        print_help_and_exit()
    
    input_tsv = os.path.expanduser(sys.argv[1])
    seq_folder = os.path.expanduser(sys.argv[2])
    output_folder = os.path.expanduser(sys.argv[3])
    use_cuda = sys.argv[4].lower() == 'true'
    compressed_matrix = sys.argv[5].lower() == 'true'
    
    # Check existence of files/folders (os.path.expanduser(sys.argv[2]))
    if not os.path.isfile(input_tsv):
        print(f"Error: Input TSV file '{input_tsv}' not found.")
        exit(1)
    
    if not os.path.isdir(seq_folder):
        print(f"Error: Sequence folder '{seq_folder}' not found.")
        exit(1)
    
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    # Checking mandatories columns in input.tsv
    required_columns = {'ID', 'model', 'scale', 'wpos', 'mpos', 'model_HFF', 'model_ESC'}
    input_tsv = pd.read_csv(input_tsv, sep='\t')
    missing_columns = required_columns - set(input_tsv.columns)
    if missing_columns:
        print(f"Error: Missing required columns in input TSV: {', '.join(missing_columns)}")
        exit(1)
        
    # Check whether columns ('model_HFF', 'model_ESC') contain only 1 or 0 :
    for column in ['model_HFF', 'model_ESC']:
        invalid_rows = input_tsv[~input_tsv[column].isin([0, 1])]
        if not invalid_rows.empty:
            invalid_ids = invalid_rows['ID'].tolist()
            print(f"Error: Column '{column}' contains values other than 0 or 1.")
            print(f"IDs concerned : {', '.join(map(str, invalid_ids))}")
            exit(1)
    
    # Check whether column ('model') contain only 32000000 or 1000000 :
    for column in ['model']:
        invalid_rows = input_tsv[~input_tsv[column].isin([32000000, 1000000])]
        if not invalid_rows.empty:
            invalid_ids = invalid_rows['ID'].tolist()
            print(f"Error: Column '{column}' contains values other than 32000000 or 1000000.")
            print(f"IDs concerned : {', '.join(map(str, invalid_ids))}")
            exit(1)
    
    #sort imput_tsv by orca model (to avoid reloading the same models)
    input_tsv = input_tsv.sort_values(by='model', ascending=False).reset_index(drop=True)
    return input_tsv, seq_folder, output_folder, use_cuda, compressed_matrix

#print + reset timer
start_time = time.time()
def print_and_reset_timer(message):
    """
    print elapsed time between each function call
    message = beginning of sentence
    """
    global start_time
    # end time
    elapsed_time = time.time() - start_time
    # Conversion to hours, minutes and seconds
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"{message} {int(hours)}h {int(minutes)}m {int(seconds)}s")
    # reset timer
    start_time = time.time()

model_dict = {32000000: '32M', 1000000: '1M'}
def load_orca_model(model, use_cuda):
    """Load ORCA model and return OrcaModel (i.e. '32M' or '1M')."""
    OrcaModel=model_dict[model]
    print(f"# Loading {OrcaModel} model...")
    orca_predict.load_resources(models=OrcaModel, use_cuda=use_cuda)
    print_and_reset_timer("## Model loaded successfully in")
    return(OrcaModel)

fasta_extensions = ['.fa', '.fa.gz', '.fa.gzip'] # list possible extentions
def read_sequence(ID, seq_dir):
    """Reads and returns the DNA sequence."""
    
    # build prefix
    fasta_prefix = os.path.join(os.path.expanduser(seq_dir), ID)
        
    # find extention
    for ext in fasta_extensions:
        fasta_path = fasta_prefix + ext
        if os.path.isfile(fasta_path):
            break
    else:
        print(f"## Error: No sequences found in '{fasta_prefix}'.")
        exit(1)
    
    # read compressed fasta
    if fasta_path.endswith(('.fa.gz', '.fa.gzip')):
        with gzip.open(fasta_path, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                print(f"## Error: No sequences found in '{fasta_path}'.")
                exit(1)
            return str(records[0].seq)
    
    elif fasta_path.endswith('.fa'):
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            print(f"## Error: No sequences found in '{fasta_path}'.")
            exit(1)
        return str(records[0].seq)
    
    else:
        print("## Error: Unsupported file format.")
        exit(1)

scale_dict = {'32000000': 0, '16000000': 1, '8000000': 2, '4000000': 3, '2000000': 4, '1000000': 5}
def export32M_matrix(row, OrcaModel, orcaCellModel, outputs, compressed_matrix):
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
      print(f"## {row['ID']} Error: Invalid scale '{e.args[0]}' in the input TSV file.")
      exit(1)
      
    param_dict = {'ID': [], 'start': [], 'end': [], 'bin_width': [], 'mpos': [], 'wpos': [], 'model': []}
    for index in scale_indexes:
      #Exports prediction and normat matrices
      matrix_pred = outputs['predictions'][0][index]
      matrix_norm = outputs['normmats'][0][index]
      suffix = list(scale_dict.keys())[index] + '_' + OrcaModel + '_' + orcaCellModel + '.tsv'

      # Save the matrices as .tsv files
      pred_file_path = os.path.join(output_dir, f"{row['ID']}_predictions_{suffix}")
      norm_file_path = os.path.join(output_dir, f"{row['ID']}_normmats_{suffix}")
      
      if compressed_matrix == True:
        with gzip.open(f"{pred_file_path}.gz", 'wt') as f_out:
            np.savetxt(f_out, matrix_pred, delimiter='\t')
        with gzip.open(f"{norm_file_path}.gz", 'wt') as f_out:
            np.savetxt(f_out, matrix_norm, delimiter='\t')
      else:
        np.savetxt(pred_file_path, matrix_pred, delimiter='\t')
        np.savetxt(norm_file_path, matrix_norm, delimiter='\t')
        
    print(f"## {OrcaModel}_{orcaCellModel} predictions exported successfully for {row['ID']}.")
    return(scale_indexes)

def export32M_parameters(row, OrcaModel, outputs, scale_indexes):
    """
    From OrcaModel 32M predictions outputs:
    export matrix parameters (start, stop...) for scale matrices to keep (as .bed).
    """
        
    param_dict = {'ID': [], 'start': [], 'end': [], 'bin_width': [], 'mpos': [], 'wpos': [], 'model': []}
    for index in scale_indexes:
      param_dict['ID'].append(row['ID'] + '_predictions_' + list(scale_dict.keys())[index])
      param_dict['start'].append(outputs['start_coords'][index])
      param_dict['end'].append(outputs['end_coords'][index])
      param_dict['bin_width'].append(str((outputs['end_coords'][index] - outputs['start_coords'][index]) / 250))
      param_dict['mpos'].append(row['mpos'])
      param_dict['wpos'].append(row['wpos'])
      param_dict['model'].append(row['model'])
      path_bed = os.path.join(output_dir, f"{row['ID']}.bed")
      pd.DataFrame.from_dict(param_dict).to_csv(path_bed, header=True, sep='\t', index=False)
    print(f"## {OrcaModel} parameters exported successfully for {row['ID']}.")

def export1M_matrix_parameters(row, OrcaModel, orcaCellModel, outputs, compressed_matrix):
    """
    From OrcaModel 1M predictions outputs:
    export matrix and parameters (start, stop...) as .bed
    """  

    # check scales values
    if not str(row['scale']) == '1000000':
      print(f"## {row['ID']} Error: Invalid scale '{row['scale']}' in the input TSV file of {row['ID']}.")
      exit(1)
    
    #Exports prediction matrices
    suffix = row['ID'] + '_predictions_' + str(row['scale']) + '_' + OrcaModel + '_' + orcaCellModel + '.tsv'
    pred_file_path = os.path.join(output_dir, f"{suffix}")

    if compressed_matrix == True:
        with gzip.open(f"{pred_file_path}.gz", 'wt') as f_out:
            np.savetxt(f_out, outputs, delimiter='\t')
    else:
        np.savetxt(pred_file_path, outputs, delimiter='\t')
    
    #Exports matrices parameters
    param_dict = {
      'ID': [row['ID'] + '_predictions_1000000'],
      'start': [row['wpos'] - 500000],
      'end': [row['wpos'] + 500000],
      'bin_width': [4000.0],
      'mpos': [row['mpos']],
      'wpos': [row['wpos']],
      'model': [row['model']]
      }
    path_bed = os.path.join(output_dir, f"{row['ID']}.bed")
    pd.DataFrame.from_dict(param_dict).to_csv(path_bed, header=True, sep='\t', index=False)
    print(f"## {OrcaModel}_{orcaCellModel} parameters exported successfully for {row['ID']}.")

def process_predictions(predictions, seq_dir, output_dir, use_cuda, compressed_matrix):
    """Processes predictions for each row of the input TSV."""
    
    #load first OrcaModel
    OrcaModel = load_orca_model(input_tsv['model'][0], use_cuda)
    
    for i, row in input_tsv.iterrows():
        # Change OrcaModel if needed
        if OrcaModel != model_dict[row['model']]:
            #end timer for the first OrcaModel prediction batch
            print_and_reset_timer(f"# {OrcaModel} ORCA predictions successfully completed in")
            #change and load new OrcaModel
            OrcaModel = load_orca_model(row['model'], use_cuda)
        
        print(f"# {row['ID']} processing...")
        
        # Validate sequence length
        sequence = read_sequence(row['ID'], seq_dir)
        if len(sequence) != row['model']:
            print(f"## Error: Sequence length does not match model size for {row['ID']}.")
            exit(1)
        
        encoded_sequence = Genome.sequence_to_encoding(sequence)[None, :, :]
        ################## 32M #########################
        if OrcaModel == '32M':
            # ORCA prediction hff
            if row['model_HFF'] == 1:
                outputs = orca_predict.genomepredict(
                    mchr=row['ID'],
                    mpos=row['mpos'],
                    wpos=row['wpos'],
                    models=['hff'],
                    sequence=encoded_sequence,
                    use_cuda=use_cuda
                    )
                # Export matrices hff
                scale_indexes = export32M_matrix(row, OrcaModel, 'hff', outputs, compressed_matrix)
            
            # ORCA prediction h1esc
            if row['model_ESC'] == 1:
                outputs = orca_predict.genomepredict(
                    mchr=row['ID'],
                    mpos=row['mpos'],
                    wpos=row['wpos'],
                    models=['h1esc'],
                    sequence=encoded_sequence,
                    use_cuda=use_cuda
                    )
                # Export matrices h1esc
                scale_indexes = export32M_matrix(row, OrcaModel, 'esc', outputs, compressed_matrix)
            
            #export parameters (same parameters for esc or hff)
            export32M_parameters(row, OrcaModel, outputs, scale_indexes)
        ###########################################
        
        ###########################################
        if OrcaModel == '1M':
            
            if use_cuda == False:
                if row['model_HFF'] == 1:
                    outputs = orca_predict.hff_1m(orca_predict.torch.FloatTensor(encoded_sequence).transpose(1, 2)).squeeze().detach().cpu().numpy()
                    export1M_matrix_parameters(row, OrcaModel, 'hff', outputs, compressed_matrix)
                
                if row['model_ESC'] == 1:
                    outputs = orca_predict.h1esc_1m(orca_predict.torch.FloatTensor(encoded_sequence).transpose(1, 2)).squeeze().detach().cpu().numpy()
                    export1M_matrix_parameters(row, OrcaModel, 'esc', outputs, compressed_matrix)
            
            if use_cuda == True:
                if row['model_HFF'] == 1:
                    outputs = orca_predict.hff_1m(orca_predict.torch.FloatTensor(encoded_sequence).cuda().transpose(1, 2)).squeeze().detach().cpu().numpy()
                    export1M_matrix_parameters(row, OrcaModel, 'esc', outputs, compressed_matrix)
                
                if row['model_ESC'] == 1:
                    outputs = orca_predict.h1esc_1m(orca_predict.torch.FloatTensor(encoded_sequence).cuda().transpose(1, 2)).squeeze().detach().cpu().numpy()
                    export1M_matrix_parameters(row, OrcaModel, 'esc', outputs, compressed_matrix)
        
    print_and_reset_timer(f"# {OrcaModel} ORCA predictions successfully completed in")


############################################
if __name__ == "__main__":
    input_tsv, seq_dir, output_dir, use_cuda, compressed_matrix = validate_inputs()
    process_predictions(input_tsv, seq_dir, output_dir, use_cuda, compressed_matrix)
