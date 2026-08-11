#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Script to predict 3D chromatin contact maps from DNA sequence using AlphaGenome.
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

# Mock imports for standalone script structure
try:
    from alphagenome.models import dna_client
    from alphagenome.data import genome
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
    parser = argparse.ArgumentParser(
        description="Predict 3D chromatin contact maps from DNA sequence using AlphaGenome.",
        epilog="""
# Input tsv file format:
------------------------
Tabulated file with header, containing as many lines as there are predictions to be made.
Mandatory columns are:
    ID: unique ID which must match fasta sequence file name: E.g, if ID = 'ID1', sequence fasta file must be 'ID1.fa' (or 'ID1.fa.gz' and 'ID1.fa.gzip').
    ontology: comma-separated list of ontology term(s) (ex: EFO:0003042,CL:0011022).
    organism: human or mouse.

# Examples:
-----------
    python alphagenomePredictionsV0.01.py -i input.tsv -f seq_folder -o output_folder --api_key YOUR_PRIVATE_KEY
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--input", "-i", type=str, required=True,
                        help="Tab-delimited file with prediction metadatas (input TSV file).")
    parser.add_argument("--fasta_dir", "-f", type=str, required=True,
                        help="Directory containing fasta sequence files.")
    parser.add_argument("--output_dir", "-o", type=str, required=True,
                        help="Directory to save output files (will be created if required).")
    parser.add_argument("--api_key", "-k", type=str, default=None,
                        help="Private/API key for AlphaGenome. Can also be set via ALPHA_GENOME_API_KEY environment variable.")
    parser.add_argument("--compress", action="store_true",
                        help="Compress predicted matrices (.gz).")
    
    return parser.parse_args()


def validate_inputs(args):
    input_tsv = os.path.expanduser(args.input)
    seq_folder = os.path.expanduser(args.fasta_dir)
    output_folder = os.path.expanduser(args.output_dir)
    
    if not os.path.isfile(input_tsv):
        logging.error(f"Input TSV file '{input_tsv}' not found.")
        sys.exit(1)
        
    if not os.path.isdir(seq_folder):
        logging.error(f"Sequence folder '{seq_folder}' not found.")
        sys.exit(1)
    
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder, exist_ok=True)
    
    # Reading input_tsv
    try:
        df = pd.read_csv(input_tsv, sep='\t')
    except Exception as e:
        logging.error(f"Failed to read input TSV: {e}")
        sys.exit(1)

    # Checking mandatory columns in input_tsv    
    required_columns = {'ID', 'ontology', 'organism'}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        logging.error(f"Missing required columns in input TSV: {', '.join(missing_columns)}")
        sys.exit(1)
    
    # Check whether organism column contains only supported values (human or mouse)
    invalid_organisms = df[~df['organism'].astype(str).str.lower().str.strip().isin(['human', 'mouse'])]
    if not invalid_organisms.empty:
        logging.error(f"Column 'organism' contains invalid values (must be 'human' or 'mouse'). IDs: {', '.join(map(str, invalid_organisms['ID'].tolist()))}")
        sys.exit(1)

    # Validate homogeneity of organism and ontology columns for batch prediction
    unique_organisms = df['organism'].astype(str).str.lower().str.strip().unique()
    if len(unique_organisms) > 1:
        logging.error(f"Batch prediction requires all rows to have the exact same organism. Found: {', '.join(unique_organisms)}")
        sys.exit(1)

    unique_ontologies = df['ontology'].astype(str).str.strip().unique()
    if len(unique_ontologies) > 1:
        logging.error("Batch prediction requires all rows to have the exact same ontology terms.")
        sys.exit(1)
        
    return df, seq_folder, output_folder


def get_api_key(args):
    api_key = args.api_key
    if not api_key:
        api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        logging.error("AlphaGenome API key not provided. Please use --api_key or set the ALPHA_GENOME_API_KEY environment variable.")
        sys.exit(1)
    return api_key


# print + reset timer
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
                        sys.exit(1)
                    return str(records[0].seq)
            else:
                with open(fasta_path, "r") as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
                    if not records:
                        logging.error(f"No sequences found in '{fasta_path}'.")
                        sys.exit(1)
                    return str(records[0].seq)
    
    logging.error(f"No sequences found for '{fasta_prefix}' in directory.")
    sys.exit(1)


# Pre-computed indices for lower triangle of 512x512 symmetric contact matrices to optimize speed
LOWER_TRIANGLE_ROWS, LOWER_TRIANGLE_COLS = np.tril_indices(512, k=-1)

def export_matrix_fast(matrix, filepath, compressed):
    matrix = np.asarray(matrix)
    formatted = matrix.astype(object)
    
    # Replace elements strictly below the diagonal with ''
    formatted[LOWER_TRIANGLE_ROWS, LOWER_TRIANGLE_COLS] = ''
    
    df_mat = pd.DataFrame(formatted)
    if compressed:
        df_mat.to_csv(filepath + '.gz', sep='\t', header=False, index=False, compression='gzip')
    else:
        df_mat.to_csv(filepath, sep='\t', header=False, index=False)


def export_metadata(dna_model, organism, ontology_terms, output_dir, metadata_filename):
    """
    Checks output metadata, filters for OutputType.CONTACT_MAPS tracks and ontology terms,
    adds a 1-based track_num column, and exports to a TSV file.
    Returns the filtered metadata DataFrame and the number of matching tracks.
    """
    try:
        from alphagenome.models import dna_client
    except ImportError:
        logging.error("Failed to import 'alphagenome' package. Please ensure it is installed.")
        sys.exit(1)

    logging.info("Checking output metadata to verify predictions...")
    try:
        output_metadata_df = dna_model.output_metadata(organism).concatenate()
    except Exception as e:
        logging.error(f"Failed to fetch output metadata: {e}")
        sys.exit(1)

    # Filter for OutputType.CONTACT_MAPS tracks and ontology terms
    contact_maps_meta = output_metadata_df[
        output_metadata_df['output_type'].astype(str).isin(['OutputType.CONTACT_MAPS', 'CONTACT_MAPS'])
    ]
    matching_tracks = contact_maps_meta[
        contact_maps_meta['ontology_curie'].isin(ontology_terms)
    ].copy()
    
    num_predictions = len(matching_tracks)
    logging.info(f"Number of interaction matrices that will be predicted per sequence: {num_predictions}")
    
    if num_predictions == 0:
        logging.error(f"No matching CONTACT_MAPS tracks found for the requested ontology terms: {ontology_terms}. Stopping execution.")
        sys.exit(1)

    # Add 1-based line/track number to metadata
    track_nums = list(range(1, num_predictions + 1))
    matching_tracks.insert(0, 'track_num', track_nums)

    # Save metadata to TSV directly
    metadata_path = os.path.join(output_dir, metadata_filename)
    matching_tracks.to_csv(metadata_path, sep='\t', index=False)
    logging.info(f"Metadata exported successfully to {metadata_path}")

    return matching_tracks, num_predictions


def process_predictions(input_tsv_df, seq_dir, output_dir, api_key, compress, matching_tracks, num_predictions):
    """Processes batch predictions for the input TSV."""
    try:
        from alphagenome.models import dna_client
        from alphagenome.data import genome
    except ImportError:
        logging.error("Failed to import 'alphagenome' package. Please ensure it is installed.")
        sys.exit(1)
        
    logging.info("Initializing AlphaGenome DNA client...")
    dna_model = dna_client.create(api_key)
    log_and_reset_timer("AlphaGenome client initialized in")
    
    # Prepare lists for batch prediction
    sequences_list = []
    intervals_list = []
    ids_list = []
    
    expected_len = 1048576
    
    for i, row in input_tsv_df.iterrows():
        ID = str(row['ID'])
        ids_list.append(ID)
        
        # Read DNA sequence
        sequence = get_sequence(ID, seq_dir)
        if len(sequence) != expected_len:
            logging.error(f"Sequence length ({len(sequence)}) does not match required size ({expected_len}) for {ID}.")
            sys.exit(1)
            
        sequences_list.append(sequence.upper())
        
        # Build interval
        interval = genome.Interval(chromosome=ID, start=1, end=1 + expected_len)
        intervals_list.append(interval)

    if not sequences_list:
        logging.warning("No sequences to process.")
        return

    # Since they are homogeneous, we can extract organism and ontology from the first row
    first_row = input_tsv_df.iloc[0]
    ontology_field = str(first_row['ontology'])
    ontology_terms = [t.strip() for t in ontology_field.split(',') if t.strip()]
    if not ontology_terms:
        logging.error("No valid ontology terms specified.")
        sys.exit(1)
        
    org_str = str(first_row['organism']).strip().lower()
    if org_str == 'human':
        organism = dna_client.Organism.HOMO_SAPIENS
    elif org_str == 'mouse':
        organism = dna_client.Organism.MUS_MUSCULUS
    else:
        logging.error(f"Unsupported organism '{org_str}'.")
        sys.exit(1)

    logging.info(f"Calling AlphaGenome predict_sequences for {len(sequences_list)} sequences ({org_str}) with terms: {ontology_terms}")
    try:
        outputs = dna_model.predict_sequences(
            sequences=sequences_list,
            organism=organism,
            requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
            ontology_terms=ontology_terms,
            intervals=intervals_list
        )
    except Exception as e:
        logging.error(f"Error calling predict_sequences: {e}")
        sys.exit(1)

    # Process and save predictions for each sequence
    for idx, ID in enumerate(ids_list):
        output = outputs[idx]
        logging.info(f"Saving contact maps for {ID}...")
        
        contact_maps = getattr(output, 'contact_maps', None)
        if contact_maps is None:
            logging.error(f"Output for {ID} does not contain 'contact_maps'.")
            sys.exit(1)
            
        if not hasattr(contact_maps, 'values'):
            logging.error(f"Contact maps object for {ID} is missing 'values' attribute.")
            sys.exit(1)
            
        values = contact_maps.values
        
        if values.ndim != 3:
            logging.error(f"Expected 3D array for contact maps, but got {values.ndim}D array with shape {values.shape}.")
            sys.exit(1)
            
        if values.shape[-1] != num_predictions:
            logging.error(f"Number of predicted tracks ({values.shape[-1]}) does not match expected metadata count ({num_predictions}) for {ID}.")
            sys.exit(1)
            
        for j in range(values.shape[-1]):
            matrice_2d = values[:, :, j]
            
            # Fetch metadata row from the passed matching_tracks
            row_meta = matching_tracks.iloc[j]
            track_num = row_meta['track_num']
            ontology_term = row_meta.get('ontology_curie')
            if pd.isna(ontology_term) or not ontology_term:
                ontology_term = row_meta.get('name', f"track_{track_num}")
                
            filename = f"{ID}_predictions_1048576_{ontology_term}_{track_num}.tsv"
            filepath = os.path.join(output_dir, filename)
            
            export_matrix_fast(matrice_2d, filepath, compress)
            
        logging.info(f"Contact maps exported successfully for {ID}.")
        
    # Reset/log the timer only when all predictions and file savings are completed
    log_and_reset_timer("Batch predictions and exports successfully completed in")


if __name__ == "__main__":
    args = parse_args()
    input_tsv_df, seq_dir, output_dir = validate_inputs(args)
    api_key = get_api_key(args)
    
    # Extract base name of the input TSV file to use in the metadata filename
    input_base = os.path.splitext(os.path.basename(args.input))[0]
    metadata_filename = f"{input_base}_AlphagenomePredictionsMetadata.tsv"
    
    # Load AlphaGenome client temporarily to fetch metadata
    try:
        from alphagenome.models import dna_client
    except ImportError:
        logging.error("Failed to import 'alphagenome' package. Please ensure it is installed.")
        sys.exit(1)
        
    # Map organism
    first_row = input_tsv_df.iloc[0]
    org_str = str(first_row['organism']).strip().lower()
    if org_str == 'human':
        organism = dna_client.Organism.HOMO_SAPIENS
    elif org_str == 'mouse':
        organism = dna_client.Organism.MUS_MUSCULUS
    else:
        logging.error(f"Unsupported organism '{org_str}'.")
        sys.exit(1)
        
    ontology_field = str(first_row['ontology'])
    ontology_terms = [t.strip() for t in ontology_field.split(',') if t.strip()]
    
    dna_model = dna_client.create(api_key)
    matching_tracks, num_predictions = export_metadata(dna_model, organism, ontology_terms, output_dir, metadata_filename)
    
    process_predictions(input_tsv_df, seq_dir, output_dir, api_key, args.compress, matching_tracks, num_predictions)
