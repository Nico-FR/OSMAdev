#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Script to predict RNA-seq tracks from DNA sequence using AlphaGenome and export as BedGraph.
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
    ALPHAGENOME_AVAILABLE = True
except ImportError:
    dna_client = None
    genome = None
    ALPHAGENOME_AVAILABLE = False

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
        description="Predict RNA-seq tracks from DNA sequence using AlphaGenome and export as BedGraph.",
        epilog="""
# Input tsv file format:
------------------------
Tabulated file with header, containing as many lines as there are predictions to be made.
Mandatory columns are:
    ID: unique ID which must match fasta sequence file name: E.g, if ID = 'ID1', sequence fasta file must be 'ID1.fa' (or 'ID1.fa.gz' and 'ID1.fa.gzip').
    ontology: comma-separated list of ontology term(s) (ex: EFO:0003042,CL:0011022).
    organism: human or mouse.

# Versions:
-----------
    1: RNAseq prediction.
    2: prediction split per 1000 predictions to avoid ram overload 

# Examples:
-----------
    python alphagenomeRNAseqV0.01.py -i input.tsv -f seq_folder -o output_folder --api_key YOUR_PRIVATE_KEY

# computing ressources:
-----------
    2.90 s per prediction and 16G of RAM.
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
                        help="Compress predicted BedGraph files (.gz).")
    
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

def check_alphagenome_installed():
    """Checks if alphagenome is installed and exits with an error if not."""
    if not ALPHAGENOME_AVAILABLE:
        logging.error("Failed to import 'alphagenome' package. Please ensure it is installed.")
        sys.exit(1)

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


def export_bedgraph(chrom, values, filepath, compressed):
    """Exports 1D track values as a single-column file containing only the value column with 5 significant digits."""

    float_format = '%.5g'
    
    if compressed:
        # Using 'wt' (write text) mode for gzip integration with numpy
        with gzip.open(filepath + '.gz', 'wt') as f:
            np.savetxt(f, values, fmt=float_format)
    else:
        np.savetxt(filepath, values, fmt=float_format)


def get_organism_and_ontology(input_tsv_df):
    """Extracts and maps organism and ontology terms from the input DataFrame."""
    check_alphagenome_installed()
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
    if not ontology_terms:
        logging.error("No valid ontology terms specified.")
        sys.exit(1)
    return organism, ontology_terms, org_str


def initialize_client(api_key):
    """Initializes and returns the AlphaGenome DNA client."""
    logging.info("Initializing AlphaGenome DNA client...")
    dna_model = dna_client.create(api_key)
    log_and_reset_timer("AlphaGenome client initialized in")
    return dna_model


def export_metadata(dna_model, organism, ontology_terms, output_dir, metadata_filename):
    """
    Checks output metadata, filters for OutputType.RNA_SEQ tracks and ontology terms,
    adds a 1-based track_num column, and exports to a TSV file.
    Returns the filtered metadata DataFrame and the number of matching tracks.
    """
    logging.info("Checking output metadata to verify predictions...")
    try:
        output_metadata_df = dna_model.output_metadata(organism).concatenate()
    except Exception as e:
        logging.error(f"Failed to fetch output metadata: {e}")
        sys.exit(1)

    # Filter for OutputType.RNA_SEQ tracks and ontology terms
    rna_seq_meta = output_metadata_df[
        output_metadata_df['output_type'].astype(str).isin(['OutputType.RNA_SEQ', 'RNA_SEQ'])
    ]
    matching_tracks = rna_seq_meta[
        rna_seq_meta['ontology_curie'].isin(ontology_terms)
    ].copy()
    
    num_predictions = len(matching_tracks)
    logging.info(f"Number of RNA-seq tracks that will be predicted per sequence: {num_predictions}")
    
    if num_predictions == 0:
        logging.error(f"No matching RNA_SEQ tracks found for the requested ontology terms: {ontology_terms}. Stopping execution.")
        sys.exit(1)

    # Add 1-based line/track number to metadata
    track_nums = list(range(1, num_predictions + 1))
    matching_tracks.insert(0, 'track_num', track_nums)

    # Save metadata to TSV directly
    metadata_path = os.path.join(output_dir, metadata_filename)
    matching_tracks.to_csv(metadata_path, sep='\t', index=False)
    logging.info(f"Metadata exported successfully to {metadata_path}")

    return matching_tracks, num_predictions


def process_predictions(dna_model, input_tsv_df, seq_dir, output_dir, compress, matching_tracks, num_predictions, organism, ontology_terms, org_str):
    """Processes batch predictions for the input TSV."""
        
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

    logging.info(f"Calling AlphaGenome predict_sequences for {len(sequences_list)} sequences ({org_str}) with terms: {ontology_terms}")
    try:
        outputs = dna_model.predict_sequences(
            sequences=sequences_list,
            organism=organism,
            requested_outputs=[dna_client.OutputType.RNA_SEQ],
            ontology_terms=ontology_terms,
            intervals=intervals_list
        )
    except Exception as e:
        logging.error(f"Error calling predict_sequences: {e}")
        sys.exit(1)

    # Process and save predictions for each sequence
    for idx, ID in enumerate(ids_list):
        output = outputs[idx]
        logging.info(f"Saving RNA-seq predictions for {ID}...")
        
        rna_seq = getattr(output, 'rna_seq', None)
        if rna_seq is None:
            logging.error(f"Output for {ID} does not contain 'rna_seq'.")
            sys.exit(1)
            
        if not hasattr(rna_seq, 'values'):
            logging.error(f"RNA-seq object for {ID} is missing 'values' attribute.")
            sys.exit(1)
            
        values = rna_seq.values
        
        if values.ndim != 2:
            logging.error(f"Expected 2D array for RNA-seq, but got {values.ndim}D array with shape {values.shape}.")
            sys.exit(1)
            
        if values.shape[-1] != num_predictions:
            logging.error(f"Number of predicted tracks ({values.shape[-1]}) does not match expected metadata count ({num_predictions}) for {ID}.")
            sys.exit(1)
            
        for j in range(values.shape[-1]):
            track_values = values[:, j]
            
            # Fetch metadata row from the passed matching_tracks
            row_meta = matching_tracks.iloc[j]
            track_num = row_meta['track_num']
            ontology_term = row_meta.get('ontology_curie')
            if pd.isna(ontology_term) or not ontology_term:
                ontology_term = row_meta.get('name', f"track_{track_num}")
                
            filename = f"{ID}_predictions_1048576_{ontology_term}_{track_num}.bedgraph"
            filepath = os.path.join(output_dir, filename)
            
            export_bedgraph(ID, track_values, filepath, compress)
            
        logging.info(f"RNA-seq predictions exported successfully for {ID}.")


def run_pipeline(args):
    """Orchestrates the entire RNA-seq prediction pipeline."""
    input_tsv_df, seq_dir, output_dir = validate_inputs(args)
    api_key = get_api_key(args)
    
    # Extract base name of the input TSV file to use in the metadata filename
    input_base = os.path.splitext(os.path.basename(args.input))[0]
    metadata_filename = f"{input_base}_AlphagenomePredictionsMetadata.tsv"
    
    organism, ontology_terms, org_str = get_organism_and_ontology(input_tsv_df)
    dna_model = initialize_client(api_key)
    
    matching_tracks, num_predictions = export_metadata(dna_model, organism, ontology_terms, output_dir, metadata_filename)
    
    # Process predictions in batches of 1000
    batch_size = 1000
    num_rows = len(input_tsv_df)

    if num_rows > batch_size:
        logging.info(f"Total number of predictions ({num_rows}) exceeds {batch_size}. Analyses will be performed in batches of {batch_size}.")
    
    for start_idx in range(0, num_rows, batch_size):
        chunk_df = input_tsv_df.iloc[start_idx : start_idx + batch_size]
        process_predictions(
            dna_model, chunk_df, seq_dir, output_dir, args.compress,
            matching_tracks, num_predictions, organism, ontology_terms, org_str
        )
        
    # Reset/log the timer only when all predictions and file savings are completed
    log_and_reset_timer("Batch predictions and exports successfully completed in")


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
