#!/usr/bin/env python3

## This script is used to:
## For each variant, count the number of isolates carrying the variant that have pDST results classified as Resistant (R) or Sensitive (S).

import pandas as pd
import argparse
import os
import sys
import logging
from collections import defaultdict

def setup_logging(output_dir, debug=False):
    """Sets up logging to file and console."""
    log_file = os.path.join(output_dir, "variant_impact_counter.log")
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Logging is set up.")

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Counts occurrences of specified variants (SNPs and Indels) in Resistant (R) and Sensitive (S) isolates "
            "based on SNP and Indel files associated with each ID."
        )
    )
    parser.add_argument('--variant_list', required=True, help='Path to the variant list file.')
    parser.add_argument('--id_list', required=True, help='Path to the ID list file.')
    parser.add_argument('--snp_dir', required=True, help='Primary directory containing SNP files named as {ID}.snp.')
    parser.add_argument('--additional_snp_dir', default=None, help='Additional directory containing SNP files named as {ID}.snp.')
    parser.add_argument('--indel_dir', required=True, help='Primary directory containing Indel files named as {ID}.indel.ano.')
    parser.add_argument('--additional_indel_dir', default=None, help='Additional directory containing Indel files named as {ID}.indel.ann.')
    parser.add_argument('--output_file', default='variant_counts.tsv', help='Path to the output TSV file.')
    parser.add_argument('--debug', action='store_true', help='Enable debug-level logging.')
    return parser.parse_args()

def load_variant_list(variant_list_path):
    """Loads the variant list into a DataFrame."""
    try:
        variant_df = pd.read_csv(
            variant_list_path,
            sep='\t',
            header=None,
            names=['variant', 'position', 'ref', 'alt', 'grading'],
            dtype={'variant': str, 'position': str, 'ref': str, 'alt': str, 'grading': int},
            engine='python'
        )
        variant_df = variant_df.dropna(subset=['variant', 'position', 'ref', 'alt'])
        # Standardize ref and alt to uppercase
        variant_df['ref'] = variant_df['ref'].str.upper().str.strip()
        variant_df['alt'] = variant_df['alt'].str.upper().str.strip()
        variant_df['position'] = variant_df['position'].str.strip()
        
        # Check for duplicate variant names with different (position, ref, alt)
        duplicate_variants = variant_df.duplicated(subset=['variant'], keep=False)
        if duplicate_variants.any():
            duplicates = variant_df[duplicate_variants]
            logging.warning(f"Found {duplicates.shape[0]} duplicate variant names with different (position, ref, alt). These will be treated as separate entries.")
        
        logging.info(f"Loaded {len(variant_df)} variants from '{variant_list_path}'.")
        return variant_df
    except Exception as e:
        logging.error(f"Error reading variant list file: {e}")
        sys.exit(1)

def load_id_list(id_list_path):
    """Loads the ID list into a DataFrame, excluding IDs with duplicated inconsistent R/S records."""
    try:
        # Read the ID list into a DataFrame
        id_df = pd.read_csv(
            id_list_path,
            sep='\t',
            header=None,
            names=['ID', 'Phenotype'],
            dtype={'ID': str, 'Phenotype': str},
            engine='python'
        )
        
        # Remove any rows with missing ID or Phenotype
        id_df.dropna(subset=['ID', 'Phenotype'], inplace=True)
        
        # Standardize the Phenotype values (ensure uppercase and no surrounding whitespace)
        id_df['Phenotype'] = id_df['Phenotype'].str.strip().str.upper()
        
        # Define valid phenotypes
        valid_pheno = {'R', 'S'}
        
        # Check for invalid phenotypes
        if not set(id_df['Phenotype']).issubset(valid_pheno):
            invalid = set(id_df['Phenotype']) - valid_pheno
            logging.error(f"Invalid phenotypes found: {invalid}. Only 'R' or 'S' are allowed.")
            sys.exit(1)
        
        # Identify duplicated IDs
        duplicated_ids = id_df[id_df.duplicated(subset=['ID'], keep=False)]
        
        # For duplicated IDs, check if all phenotypes are the same
        phenotype_counts = duplicated_ids.groupby('ID')['Phenotype'].nunique()
        inconsistent_ids = phenotype_counts[phenotype_counts > 1].index.tolist()
        
        # Log and exclude inconsistent IDs
        if inconsistent_ids:
            id_df = id_df[~id_df['ID'].isin(inconsistent_ids)]
            logging.warning(f"Excluded {len(inconsistent_ids)} IDs due to inconsistent Phenotype: {inconsistent_ids}")
        else:
            logging.info("No inconsistent Phenotypes found among duplicated IDs.")
        
        # Remove duplicated IDs, keeping the first occurrence
        id_df = id_df.drop_duplicates(subset=['ID'], keep='first')
        
        logging.info(f"Loaded {len(id_df)} unique IDs with consistent Phenotypes.")
        return id_df
    except Exception as e:
        logging.error(f"Error reading ID list file: {e}")
        sys.exit(1)

def preload_snp_files(id_df, snp_dir, additional_snp_dir=None):
    """
    Preloads all SNP files into a dictionary mapping (position, ref, alt) to set of IDs.
    Returns the dictionary and the set of IDs that have SNP files.
    """
    snp_variants = defaultdict(set)
    ids_with_snp = set()
    missing_snp = 0
    total_ids = len(id_df)
    
    for idx, row in id_df.iterrows():
        id = row['ID']
        snp_file_primary = os.path.join(snp_dir, f"{id}.snp")
        snp_file_additional = os.path.join(additional_snp_dir, f"{id}.snp") if additional_snp_dir else None
        
        if os.path.isfile(snp_file_primary):
            snp_file = snp_file_primary
        elif snp_file_additional and os.path.isfile(snp_file_additional):
            snp_file = snp_file_additional
        else:
            missing_snp += 1
            continue
        
        try:
            snp_df = pd.read_csv(
                snp_file,
                sep='\t',
                header=None,
                names=['position', 'ref', 'alt'],
                dtype={'position': str, 'ref': str, 'alt': str},
                engine='python'
            )
            snp_df = snp_df.dropna(subset=['position', 'ref', 'alt'])
            snp_df['position'] = snp_df['position'].str.strip()
            snp_df['ref'] = snp_df['ref'].str.upper().str.strip()
            snp_df['alt'] = snp_df['alt'].str.upper().str.strip()
            
            for _, snp_row in snp_df.iterrows():
                variant_key = (snp_row['position'], snp_row['ref'], snp_row['alt'])
                snp_variants[variant_key].add(id)
            
            ids_with_snp.add(id)
        except Exception as e:
            logging.warning(f"Error reading SNP file '{snp_file}': {e}.")
            missing_snp += 1
            continue
        
        # Optional: Progress logging every 1000 IDs
        if (idx + 1) % 1000 == 0 or (idx + 1) == total_ids:
            logging.info(f"Preloaded {idx + 1}/{total_ids} SNP files.")
    
    logging.info(f"Preloaded SNP data for {len(snp_variants)} unique SNP variants.")
    logging.info(f"Total IDs with SNP files: {len(ids_with_snp)}")
    logging.info(f"Total IDs excluded due to missing or unreadable SNP files: {missing_snp}.")
    return snp_variants, ids_with_snp, missing_snp

def preload_indel_files(id_df, indel_dir, additional_indel_dir=None):
    """
    Preloads all Indel files into a dictionary mapping (position, ref, alt) to set of IDs.
    Returns the dictionary and the set of IDs that have Indel files.
    """
    indel_variants = defaultdict(set)
    ids_with_indel = set()
    missing_indel = 0
    total_ids = len(id_df)
    
    for idx, row in id_df.iterrows():
        id = row['ID']
        indel_file_primary = os.path.join(indel_dir, f"{id}.indel.ano")
        indel_file_additional = os.path.join(additional_indel_dir, f"{id}.indel.ann") if additional_indel_dir else None
        
        if os.path.isfile(indel_file_primary):
            indel_file = indel_file_primary
        elif indel_file_additional and os.path.isfile(indel_file_additional):
            indel_file = indel_file_additional
        else:
            missing_indel += 1
            continue
        
        try:
            indel_df = pd.read_csv(
                indel_file,
                sep='\t',
                header=None,
                usecols=[0, 1, 2],  # Only first three columns
                names=['position', 'ref', 'alt'],
                dtype={'position': str, 'ref': str, 'alt': str},
                engine='python'
            )
            indel_df = indel_df.dropna(subset=['position', 'ref', 'alt'])
            indel_df['position'] = indel_df['position'].str.strip()
            indel_df['ref'] = indel_df['ref'].str.upper().str.strip()
            indel_df['alt'] = indel_df['alt'].str.upper().str.strip()
            
            for _, indel_row in indel_df.iterrows():
                variant_key = (indel_row['position'], indel_row['ref'], indel_row['alt'])
                indel_variants[variant_key].add(id)
            
            ids_with_indel.add(id)
        except Exception as e:
            logging.warning(f"Error reading Indel file '{indel_file}': {e}.")
            missing_indel += 1
            continue
        
        # Optional: Progress logging every 1000 IDs
        if (idx + 1) % 1000 == 0 or (idx + 1) == total_ids:
            logging.info(f"Preloaded {idx + 1}/{total_ids} Indel files.")
    
    logging.info(f"Preloaded Indel data for {len(indel_variants)} unique Indel variants.")
    logging.info(f"Total IDs with Indel files: {len(ids_with_indel)}")
    logging.info(f"Total IDs excluded due to missing or unreadable Indel files: {missing_indel}.")
    return indel_variants, ids_with_indel, missing_indel

def count_variants(variant_df, snp_variants, indel_variants, ids_final, id_df):
    """
    Counts occurrences of each variant in R and S isolates.
    Each variant is uniquely identified by (variant name, position, ref, alt).
    """
    # Initialize count columns
    variant_df['R_count'] = 0
    variant_df['S_count'] = 0
    
    # Create a mapping from (position, ref, alt) to set of IDs (from both SNP and Indel)
    combined_variants = defaultdict(set)
    
    for variant_key, ids in snp_variants.items():
        combined_variants[variant_key].update(ids)
    
    for variant_key, ids in indel_variants.items():
        combined_variants[variant_key].update(ids)
    
    logging.info(f"Total unique combined variants (SNPs + Indels): {len(combined_variants)}")
    
    # Iterate through each variant row and assign counts
    for idx, row in variant_df.iterrows():
        variant_name = row['variant']
        position = row['position']
        ref = row['ref']
        alt = row['alt']
        variant_key = (position, ref, alt)
        
        ids_with_variant = combined_variants.get(variant_key, set()).intersection(ids_final)
        
        if not ids_with_variant:
            # No IDs have this variant
            variant_df.at[idx, 'R_count'] = 0
            variant_df.at[idx, 'S_count'] = 0
            continue
        
        # Fetch phenotypes for these IDs
        phenotypes = id_df.set_index('ID').loc[list(ids_with_variant), 'Phenotype']
        
        # Count R and S
        r_count = (phenotypes == 'R').sum()
        s_count = (phenotypes == 'S').sum()
        
        # Assign counts
        variant_df.at[idx, 'R_count'] = r_count
        variant_df.at[idx, 'S_count'] = s_count
        
        # Optional: Progress logging every 1000 variants
        if (idx + 1) % 1000 == 0 or (idx + 1) == len(variant_df):
            logging.info(f"Processed {idx + 1}/{len(variant_df)} variants.")
    
    logging.info("Completed counting variants.")
    return variant_df

def save_results(variant_df, output_file):
    """Saves the variant counts to a TSV file."""
    # Save as tab-separated file
    variant_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Saved variant counts to '{output_file}' as a tab-separated file.")

def main():
    args = parse_arguments()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Setup logging with debug flag
    setup_logging(output_dir if output_dir else '.', debug=args.debug)
    
    # Load data
    variant_df = load_variant_list(args.variant_list)
    id_df = load_id_list(args.id_list)
    
    # Preload SNP and Indel files
    snp_variants, ids_with_snp, missing_snp = preload_snp_files(id_df, args.snp_dir, args.additional_snp_dir)
    indel_variants, ids_with_indel, missing_indel = preload_indel_files(id_df, args.indel_dir, args.additional_indel_dir)
    
    # Find IDs that have both SNP and Indel files
    ids_final = ids_with_snp.intersection(ids_with_indel)
    excluded_ids = len(id_df) - len(ids_final)
    logging.info(f"Total IDs included (have both SNP and Indel files): {len(ids_final)}")
    logging.info(f"Total IDs excluded (do not have both SNP and Indel files): {excluded_ids}")
    
    if not ids_final:
        logging.error("No IDs have both SNP and Indel files. Exiting.")
        sys.exit(1)
    
    # Count variants
    variant_df = count_variants(variant_df, snp_variants, indel_variants, ids_final, id_df)
    
    # Save results
    save_results(variant_df, args.output_file)
    
    # Report missing files
    total_ids = len(id_df)
    included_ids = len(ids_final)
    logging.info(f"Total IDs processed: {included_ids}")
    logging.info(f"Total IDs excluded due to missing SNP files: {missing_snp}")
    logging.info(f"Total IDs excluded due to missing Indel files: {missing_indel}")
    
    print(f"Total IDs processed: {included_ids}")
    print(f"Total IDs excluded due to missing SNP files: {missing_snp}")
    print(f"Total IDs excluded due to missing Indel files: {missing_indel}")

if __name__ == "__main__":
    main()
