#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Script to predict 3D organization from DNA sequence using AlphaGenome.
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

# Mock imports for AlphaGenome (so the script can be written and parsed)
try:
    from alphagenome import dna_client
    from alphagenome import dna_model
    from alphagenome.genome import Interval
except ImportError:
    pass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

def parse_args():
    epilog="""
# Input tsv file format:
------------------------
Tabulated file with header, containing as many lines as there are predictions to be made.
Mandatory columns are:
    ID: unique ID which must match fasta sequence file: E.g, if ID = 'ID1', sequence fasta file must be 'ID1.fa'.
    model: 1000000 (AlphaGenome handles 1MB inputs natively).
    scale: 1048576.
    start.mat: The start coordinate of the 1048576 bp sequence.
    stop.mat: The stop coordinate of the 1048576 bp sequence.
    model_HFF: 0 or 1 to save predictions from HFF model ('EFO:0001196' - Human foreskin fibroblast).
    model_ESC: 0 or 1 to save predictions from ESC model ('EFO:0003042' - H1 human embryonic stem cell line).

# Examples:
-----------
    python alphaGenomePredictions.py -i input.tsv -f seq_folder -o output_folder --compress
    """
    parser = argparse.ArgumentParser(description="Predict 3D organization from DNA sequence using AlphaGenome.",
                                     epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--input", "-i", type=str, required=True,
                        help="Tab-delimited file with prediction parameters (input TSV file).")
    parser.add_argument("--fasta_dir", "-f", type=str, required=True,
                        help="Directory containing fasta sequence files.")
    parser.add_argument("--output", "-o", type=str, required=True,
                        help="Directory to save output files (will be created if required).")
    parser.add_argument("--compress", action="store_true",
                        help="Compress predicted matrices (.gz).")

    return parser.parse_args()


def validate_inputs(args):
    input_tsv = os.path.expanduser(args.input)
    seq_folder = os.path.expanduser(args.fasta_dir)
    output_folder = os.path.expanduser(args.output)

    if not os.path.isfile(input_tsv):
        logging.error(f"Input TSV file '{input_tsv}' not found.")
        exit(1)

    if not os.path.isdir(seq_folder):
        logging.error(f"Sequence folder '{seq_folder}' not found.")
        exit(1)

    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    required_columns = {'ID', 'model', 'scale', 'start.mat', 'stop.mat', 'model_HFF', 'model_ESC'}
    try:
        df = pd.read_csv(input_tsv, sep='\t')
    except Exception as e:
        logging.error(f"Failed to read input TSV: {e}")
        exit(1)

    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        logging.error(f"Missing required columns in input TSV: {', '.join(missing_columns)}")
        exit(1)

    for column in ['model_HFF', 'model_ESC']:
        invalid_rows = df[~df[column].isin([0, 1])]
        if not invalid_rows.empty:
            logging.error(f"Column '{column}' contains values other than 0 or 1. IDs: {', '.join(map(str, invalid_rows['ID'].tolist()))}")
            exit(1)

    # Ensure all models are 1MB
    invalid_rows = df[~df['model'].isin([1048576, 1000000])]
    if not invalid_rows.empty:
        logging.warning(f"AlphaGenome expects 1MB (1000000) model lengths. Found other values for IDs: {', '.join(map(str, invalid_rows['ID'].tolist()))}. Forcing 1MB.")
        df['model'] = 1048576

    return df, seq_folder, output_folder


start_time = time.time()
def log_and_reset_timer(message):
    global start_time
    elapsed_time = time.time() - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    logging.info(f"{message} {int(hours)}h {int(minutes)}m {int(seconds)}s")
    start_time = time.time()


fasta_extensions = ['.fa', '.fa.gz', '.fa.gzip']
def get_sequence(ID, seq_dir):
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


def export_matrix_fast(matrix, filepath, compressed):
    df_mat = pd.DataFrame(matrix)
    if compressed:
        df_mat.to_csv(filepath + '.gz', sep='\t', header=False, index=False, compression='gzip')
    else:
        df_mat.to_csv(filepath, sep='\t', header=False, index=False)


def export1M_matrix_parameters(row, CellModelStr, matrix, compressed_matrix, output_dir):
    suffix = str(row['ID']) + f'_predictions_1048576_{CellModelStr}.tsv'
    pred_file_path = os.path.join(output_dir, suffix)

    export_matrix_fast(matrix, pred_file_path, compressed_matrix)

    param_dict = {
      'ID': [str(row['ID']) + '_predictions_1048576'],
      'start': [row['start.mat']],
      'end': [row['stop.mat']],
      'bin_width': [2048.0],
      'model': [row['model']]
    }
    path_bed = os.path.join(output_dir, f"{row['ID']}.bed")
    pd.DataFrame.from_dict(param_dict).to_csv(path_bed, header=True, sep='\t', index=False)
    logging.info(f"AlphaGenome_{CellModelStr} parameters exported successfully for {row['ID']}.")


def process_predictions(input_tsv_df, seq_dir, output_dir, compressed_matrix):
    # AlphaGenome mapping for required ontologies
    # Orca HFF => EFO:0001196 (Human foreskin fibroblast)
    # Orca ESC => EFO:0003042 (H1 human embryonic stem cell line)

    for i, row in input_tsv_df.iterrows():
        logging.info(f"Processing {row['ID']}...")

        sequence_full = get_sequence(row['ID'], seq_dir)

        # AlphaGenome specifically takes exactly 1048576 bp
        target_len = 1048576
        seq_len = len(sequence_full)

        if seq_len < target_len:
            logging.error(f"Sequence length ({seq_len}) is smaller than AlphaGenome's required size ({target_len}) for {row['ID']}.")
            exit(1)

        # If sequence is larger, we take the center piece to mimic Orca behavior
        start_idx = (seq_len - target_len) // 2
        end_idx = start_idx + target_len
        sequence = sequence_full[start_idx:end_idx].upper()

        start_pos = int(row['start.mat'])
        end_pos = int(row['stop.mat'])

        # Create Interval
        interval = Interval(chromosome=str(row['ID']), start=start_pos, end=end_pos)

        ontologies = []
        if row['model_HFF'] == 1:
            ontologies.append('EFO:0001196')
        if row['model_ESC'] == 1:
            ontologies.append('EFO:0003042')

        if not ontologies:
            logging.warning(f"No cell model requested for {row['ID']}. Skipping prediction.")
            continue

        logging.info(f"Predicting contact maps for {row['ID']} using ontologies {ontologies}...")

        output = dna_model.predict_sequence(
            sequence,
            requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
            ontology_terms=ontologies,
            interval=interval
        )

        contact_maps = output.contact_maps

        if not hasattr(contact_maps, 'values') or not hasattr(contact_maps, 'metadata'):
            logging.error(f"Prediction output for {row['ID']} does not contain values or metadata.")
            continue

        values = contact_maps.values
        metadata_df = contact_maps.metadata

        # Determine track names
        try:
            track_names = metadata_df['name'].tolist()
        except KeyError:
            logging.warning(f"Column 'name' not found in metadata for {row['ID']}. Using fallback.")
            track_names = [f"track_{idx}" for idx in range(values.shape[-1])]

        # Verify 3D dimensions
        if values.ndim != 3:
            logging.error(f"Contact maps expected 3D array, got {values.ndim}D for {row['ID']}.")
            continue

        # Export corresponding matrices
        for idx, name in enumerate(track_names):
            matrice_2d = values[:, :, idx]

            # Identify which ontology this track matches to route to correct ESC/HFF naming
            cell_str = "unknown"
            if "EFO:0001196" in name or "HFF" in name.upper() or "fibroblast" in name.lower():
                cell_str = "hff"
            elif "EFO:0003042" in name or "ESC" in name.upper() or "embryonic" in name.lower():
                cell_str = "esc"
            else:
                # Fallback to order if we can't match string exactly, though metadata should have the term
                cell_str = "hff" if row['model_HFF'] == 1 and idx == 0 else "esc"

            export1M_matrix_parameters(row, cell_str, matrice_2d, compressed_matrix, output_dir)

    log_and_reset_timer("AlphaGenome predictions successfully completed in")


if __name__ == "__main__":
    args = parse_args()
    input_tsv_df, seq_dir, output_dir = validate_inputs(args)
    # Note: AlphaGenome doesn't seem to explicitly require a --cuda switch in its python API
    # predict_sequence wrapper, it likely handles the device dynamically under the hood, so
    # we omitted args.cuda from the predict flow, though it can be added back if needed.
    process_predictions(input_tsv_df, seq_dir, output_dir, args.compress)
