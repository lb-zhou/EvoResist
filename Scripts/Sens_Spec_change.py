#!/usr/bin/env python3

## This script is used to:
## For each variant, compute the change in prediction sensitivity and specificity compared to models using mutations from WHO G1 or WHO G1 + G2.

import pandas as pd
import argparse
import os
import sys
import logging
import math
from collections import defaultdict

def setup_logging(output_dir, debug=False):
    """Sets up logging to file and console."""
    log_file = os.path.join(output_dir, "variant_impact_analysis.log")
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
            "Evaluate the impact of each variant in List2 on Sensitivity, Specificity, and Precision "
            "when added to List1, based on SNP and Indel data associated with each ID."
        )
    )
    parser.add_argument('--list1_variants_file', required=True, help='Path to the List1 variants file.')
    parser.add_argument('--list2_variants_file', required=True, help='Path to the List2 variants file.')
    parser.add_argument('--id_list_file', required=True, help='Path to the ID list file.')
    parser.add_argument('--snp_dir', required=True, help='Primary directory containing SNP files named as {ID}.snp.')
    parser.add_argument('--additional_snp_dir', default=None, help='Additional directory containing SNP files named as {ID}.snp.')
    parser.add_argument('--indel_dir', required=True, help='Primary directory containing Indel files named as {ID}.indel.ano.')
    parser.add_argument('--additional_indel_dir', default=None, help='Additional directory containing Indel files named as {ID}.indel.ann.')
    parser.add_argument('--output_dir', default='output_impact_analysis', help='Directory to save output files and logs.')
    parser.add_argument('--debug', action='store_true', help='Enable debug-level logging.')
    return parser.parse_args()

def load_list1_variants(list1_file):
    """Loads List1 variants into a DataFrame."""
    try:
        list1_df = pd.read_csv(
            list1_file,
            sep='\t',
            header=None,
            names=['variant', 'position', 'ref', 'alt', 'grading'],
            dtype={'variant': str, 'position': str, 'ref': str, 'alt': str, 'grading': int},
            engine='python'
        )
        list1_df.dropna(subset=['variant', 'position', 'ref', 'alt'], inplace=True)
        # Standardize formatting
        list1_df['ref'] = list1_df['ref'].str.upper().str.strip()
        list1_df['alt'] = list1_df['alt'].str.upper().str.strip()
        list1_df['position'] = list1_df['position'].str.strip()
        logging.info(f"Loaded {len(list1_df)} variants from List1.")
        return list1_df
    except Exception as e:
        logging.error(f"Error reading List1 variants file: {e}")
        sys.exit(1)

def load_list2_variants(list2_file):
    """Loads List2 variants into a DataFrame."""
    try:
        list2_df = pd.read_csv(
            list2_file,
            sep='\t',
            header=None,
            names=['variant', 'grading', 'position', 'ref', 'alt', 'convergent_time'],
            dtype={'variant': str, 'grading': int, 'position': str, 'ref': str, 'alt': str, 'convergent_time': int},
            engine='python'
        )
        list2_df.dropna(subset=['variant', 'position', 'ref', 'alt'], inplace=True)
        # Standardize formatting
        list2_df['ref'] = list2_df['ref'].str.upper().str.strip()
        list2_df['alt'] = list2_df['alt'].str.upper().str.strip()
        list2_df['position'] = list2_df['position'].str.strip()
        logging.info(f"Loaded {len(list2_df)} variants from List2.")
        return list2_df
    except Exception as e:
        logging.error(f"Error reading List2 variants file: {e}")
        sys.exit(1)

def load_id_list(id_list_file):
    """Loads the ID list into a DataFrame, excluding IDs with duplicated inconsistent R/S records."""
    try:
        # Read the ID list into a DataFrame
        id_df = pd.read_csv(
            id_list_file,
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

def create_variant_mapping(list1_df, list2_df):
    """
    Creates a mapping from (variant name, position, ref, alt) to a unique identifier.
    This ensures that variants with the same name but different positions or alleles are treated separately.
    """
    variant_map = {}
    # Add List1 variants
    for _, row in list1_df.iterrows():
        key = (row['variant'], row['position'], row['ref'], row['alt'])
        variant_map[key] = {
            'variant': row['variant'],
            'position': row['position'],
            'ref': row['ref'],
            'alt': row['alt']
        }
    # Add List2 variants
    for _, row in list2_df.iterrows():
        key = (row['variant'], row['position'], row['ref'], row['alt'])
        variant_map[key] = {
            'variant': row['variant'],
            'position': row['position'],
            'ref': row['ref'],
            'alt': row['alt']
        }
    logging.info("Created combined variant mapping from List1 and List2.")
    return variant_map

def load_all_snp_files(id_df, snp_dir, additional_snp_dir=None):
    """
    Loads all SNP files into a dictionary {ID: set of (position, ref, alt)}.
    Returns the dictionary and the count of missing or unreadable SNP files.
    """
    snp_dict = {}
    missing_snp = 0
    total_ids = len(id_df)
    for idx, row in id_df.iterrows():
        id = row['ID']
        # Check in primary snp_dir
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
                usecols=[0,1,2],  # === MODIFIED ===
                dtype={'position': str, 'ref': str, 'alt': str},
                engine='python'
            )
            snp_df.dropna(subset=['position', 'ref', 'alt'], inplace=True)
            if snp_df.empty:
                missing_snp += 1
                continue
            # Standardize formatting
            snp_df['position'] = snp_df['position'].str.strip()
            snp_df['ref'] = snp_df['ref'].str.upper().str.strip()
            snp_df['alt'] = snp_df['alt'].str.upper().str.strip()
            snp_set = set(zip(snp_df['position'], snp_df['ref'], snp_df['alt']))
            snp_dict[id] = snp_set
        except Exception as e:
            logging.warning(f"Error reading SNP file '{snp_file}': {e}.")
            missing_snp += 1
            continue
        # Optional: Progress logging every 1000 IDs
        if (idx + 1) % 1000 == 0 or (idx + 1) == total_ids:
            logging.info(f"Processed {idx + 1}/{total_ids} SNP files.")
    logging.info(f"Loaded SNP data for {len(snp_dict)} IDs. Missing or unreadable SNP files: {missing_snp}.")
    return snp_dict, missing_snp

def load_all_indel_files(id_df, indel_dir, additional_indel_dir=None):
    """
    Loads all Indel files into a dictionary {ID: set of (position, ref, alt)}.
    Returns the dictionary and the count of missing or unreadable Indel files.
    """
    indel_dict = {}
    missing_indel = 0
    total_ids = len(id_df)
    for idx, row in id_df.iterrows():
        id = row['ID']
        # Check in primary indel_dir
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
                names=['position', 'ref', 'alt'],
                usecols=[0,1,2],  # === MODIFIED ===
                dtype={'position': str, 'ref': str, 'alt': str},
                engine='python'
            )
            indel_df.dropna(subset=['position', 'ref', 'alt'], inplace=True)
            if indel_df.empty:
                missing_indel += 1
                continue
            # Standardize formatting
            indel_df['position'] = indel_df['position'].str.strip()
            indel_df['ref'] = indel_df['ref'].str.upper().str.strip()
            indel_df['alt'] = indel_df['alt'].str.upper().str.strip()
            indel_set = set(zip(indel_df['position'], indel_df['ref'], indel_df['alt']))
            indel_dict[id] = indel_set
        except Exception as e:
            logging.warning(f"Error reading Indel file '{indel_file}': {e}.")
            missing_indel += 1
            continue
        # Optional: Progress logging every 1000 IDs
        if (idx + 1) % 1000 == 0 or (idx + 1) == total_ids:
            logging.info(f"Processed {idx + 1}/{total_ids} Indel files.")
    logging.info(f"Loaded Indel data for {len(indel_dict)} IDs. Missing or unreadable Indel files: {missing_indel}.")
    return indel_dict, missing_indel

def generate_baseline_predictions(id_df, snp_dict, indel_dict, list1_set):
    """
    Generates baseline predictions using only List1 variants.
    Returns a DataFrame with predictions and the baseline metrics.
    """
    predictions = []
    phenotypes = []
    processed_ids = 0
    total_ids = len(snp_dict)
    for id in snp_dict.keys():  # Iterate only over IDs with SNP and Indel data
        phenotype = id_df.loc[id_df['ID'] == id, 'Phenotype'].values[0]
        snp_set = snp_dict.get(id, set())
        indel_set = indel_dict.get(id, set())
        combined_set = snp_set.union(indel_set)
        has_list1 = bool(combined_set.intersection(list1_set))
        prediction = 'R' if has_list1 else 'S'
        predictions.append(prediction)
        phenotypes.append(phenotype)
        processed_ids += 1
        if (processed_ids) % 1000 == 0 or (processed_ids) == total_ids:
            logging.info(f"Generated baseline predictions for {processed_ids}/{total_ids} IDs.")
    predictions_df = pd.DataFrame({
        'ID': list(snp_dict.keys()),
        'Phenotype': phenotypes,
        'Prediction': predictions
    })
    baseline_metrics = compute_metrics(predictions_df)
    logging.info(f"Baseline Metrics: Sensitivity={baseline_metrics['Sensitivity']}, Specificity={baseline_metrics['Specificity']}, Precision={baseline_metrics['Precision']}")
    return predictions_df, baseline_metrics

def are_close(a, b, tol=1e-6):
    """Checks if two floating-point numbers are close within a specified tolerance."""
    return math.isclose(a, b, abs_tol=tol)

def compute_metrics(df):
    """Computes Sensitivity, Specificity, and Precision based on predictions."""
    TP = ((df['Phenotype'] == 'R') & (df['Prediction'] == 'R')).sum()
    TN = ((df['Phenotype'] == 'S') & (df['Prediction'] == 'S')).sum()
    FP = ((df['Phenotype'] == 'S') & (df['Prediction'] == 'R')).sum()
    FN = ((df['Phenotype'] == 'R') & (df['Prediction'] == 'S')).sum()

    Sensitivity = TP / (TP + FN) if (TP + FN) > 0 else None
    Specificity = TN / (TN + FP) if (TN + FP) > 0 else None
    Precision = TP / (TP + FP) if (TP + FP) > 0 else None

    metrics = {
        'TP': TP,
        'TN': TN,
        'FP': FP,
        'FN': FN,
        'Sensitivity': Sensitivity,
        'Specificity': Specificity,
        'Precision': Precision
    }

    return metrics

def calculate_metrics(new_R_ids, new_S_ids, actual_R_ids, actual_S_ids):
    """Calculates Sensitivity, Specificity, and Precision based on new R and S IDs."""
    TP = len(new_R_ids.intersection(actual_R_ids))
    TN = len(new_S_ids.intersection(actual_S_ids))
    FP = len(new_R_ids.intersection(actual_S_ids))
    FN = len(new_S_ids.intersection(actual_R_ids))

    Sensitivity = TP / (TP + FN) if (TP + FN) > 0 else None
    Specificity = TN / (TN + FP) if (TN + FP) > 0 else None
    Precision = TP / (TP + FP) if (TP + FP) > 0 else None

    return Sensitivity, Specificity, Precision

def analyze_variant_impact(list2_df, list1_set, snp_dict, indel_dict, baseline_metrics, id_df, output_dir):
    """
    Analyzes the impact of each List2 variant by adding it to List1 and computing metric changes.
    Returns a DataFrame with impact metrics.
    """
    impact_results = []
    total_variants = len(list2_df)
    processed_variants = 0

    # Corrected Definitions: Restrict to included IDs only
    actual_R_ids = set(id_df[(id_df['Phenotype'] == 'R') & (id_df['ID'].isin(snp_dict.keys()))]['ID'])
    actual_S_ids = set(id_df[(id_df['Phenotype'] == 'S') & (id_df['ID'].isin(snp_dict.keys()))]['ID'])

    # Precompute IDs with List1 variants
    ids_with_list1 = set()
    for id, snp_set in snp_dict.items():
        indel_set = indel_dict.get(id, set())
        combined_set = snp_set.union(indel_set)
        if combined_set.intersection(list1_set):
            ids_with_list1.add(id)
    logging.debug(f"IDs with List1 variants: {ids_with_list1}")

    # Precompute mapping from variant to IDs that have it (both SNP and Indel)
    variant_to_ids = defaultdict(set)
    for id, snp_set in snp_dict.items():
        indel_set = indel_dict.get(id, set())
        combined_set = snp_set.union(indel_set)
        for variant in combined_set:
            variant_to_ids[variant].add(id)
    logging.debug(f"Precomputed variant to IDs mapping.")

    for _, variant_row in list2_df.iterrows():
        variant = variant_row['variant']
        position = variant_row['position']
        ref = variant_row['ref']
        alt = variant_row['alt']
        list2_variant = (variant, position, ref, alt)  # Changed to include variant name

        logging.info(f"Analyzing impact of variant '{variant}' ({position}, {ref}, {alt}).")

        # Get IDs that have this variant
        variant_present_ids = variant_to_ids.get((position, ref, alt), set())
        logging.debug(f"Variant '{variant}' present in {len(variant_present_ids)} isolates: {variant_present_ids}")

        # Identify IDs that change prediction from S to R
        # These are IDs that have the variant and were previously predicted as S (i.e., not in ids_with_list1)
        ids_changed_S_to_R = variant_present_ids - ids_with_list1
        logging.debug(f"Variant '{variant}' causes {len(ids_changed_S_to_R)} IDs to change from S to R.")

        # Determine if the variant has no presence in any isolate
        is_no_presence = not variant_present_ids
        if is_no_presence:
            # Metrics should remain unchanged
            Sensitivity = baseline_metrics['Sensitivity']
            Specificity = baseline_metrics['Specificity']
            Precision = baseline_metrics['Precision']
            delta_sensitivity = 0.0
            delta_specificity = 0.0
            delta_precision = 0.0
            logging.debug(f"Variant '{variant}' is not present in any isolates. Metrics remain unchanged.")

            # Append the results
            impact_results.append({
                'Variant': variant,
                'Position': position,
                'Ref': ref,
                'Alt': alt,
                'Sensitivity': round(Sensitivity, 6) if Sensitivity is not None else None,
                'Delta_Sensitivity': delta_sensitivity,
                'Specificity': round(Specificity, 6) if Specificity is not None else None,
                'Delta_Specificity': delta_specificity,
                'Precision': round(Precision, 6) if Precision is not None else None,
                'Delta_Precision': delta_precision,
                'IDs_changed_S_to_R': ','.join(sorted(ids_changed_S_to_R)) if ids_changed_S_to_R else ''
            })

            # Assertions with tolerance
            if Sensitivity is not None and baseline_metrics['Sensitivity'] is not None:
                assert are_close(Sensitivity, baseline_metrics['Sensitivity']), f"Sensitivity mismatch for variant '{variant}'. Expected {baseline_metrics['Sensitivity']}, got {Sensitivity}."
            if Specificity is not None and baseline_metrics['Specificity'] is not None:
                assert are_close(Specificity, baseline_metrics['Specificity']), f"Specificity mismatch for variant '{variant}'. Expected {baseline_metrics['Specificity']}, got {Specificity}."
            if Precision is not None and baseline_metrics['Precision'] is not None:
                assert are_close(Precision, baseline_metrics['Precision']), f"Precision mismatch for variant '{variant}'. Expected {baseline_metrics['Precision']}, got {Precision}."
        else:
            # Variant is present in at least one isolate (could be R or S)
            # Compute new R and S IDs
            new_R_ids = ids_with_list1.union(variant_present_ids)
            new_S_ids = set(snp_dict.keys()) - new_R_ids

            logging.debug(f"New R IDs count: {len(new_R_ids)}")
            logging.debug(f"New S IDs count: {len(new_S_ids)}")

            # Calculate TP, TN, FP, FN using set operations
            TP = len(new_R_ids.intersection(actual_R_ids))
            TN = len(new_S_ids.intersection(actual_S_ids))
            FP = len(new_R_ids.intersection(actual_S_ids))
            FN = len(new_S_ids.intersection(actual_R_ids))

            logging.debug(f"TP: {TP}, TN: {TN}, FP: {FP}, FN: {FN}")

            # Compute new metrics
            Sensitivity, Specificity, Precision = calculate_metrics(new_R_ids, new_S_ids, actual_R_ids, actual_S_ids)
            Sensitivity = round(Sensitivity, 6) if Sensitivity is not None else None
            Specificity = round(Specificity, 6) if Specificity is not None else None
            Precision = round(Precision, 6) if Precision is not None else None

            logging.debug(f"Sensitivity: {Sensitivity}, Specificity: {Specificity}, Precision: {Precision}")

            # Compute delta metrics
            delta_sensitivity = (Sensitivity - baseline_metrics['Sensitivity']) if (Sensitivity is not None and baseline_metrics['Sensitivity'] is not None) else None
            delta_specificity = (Specificity - baseline_metrics['Specificity']) if (Specificity is not None and baseline_metrics['Specificity'] is not None) else None
            delta_precision = (Precision - baseline_metrics['Precision']) if (Precision is not None and baseline_metrics['Precision'] is not None) else None

            logging.debug(f"Delta Sensitivity: {delta_sensitivity}, Delta Specificity: {delta_specificity}, Delta Precision: {delta_precision}")

            # Append the results, including the list of IDs that changed from S to R
            impact_results.append({
                'Variant': variant,
                'Position': position,
                'Ref': ref,
                'Alt': alt,
                'Sensitivity': Sensitivity,
                'Delta_Sensitivity': delta_sensitivity,
                'Specificity': Specificity,
                'Delta_Specificity': delta_specificity,
                'Precision': Precision,
                'Delta_Precision': delta_precision,
                'IDs_changed_S_to_R': ','.join(sorted(ids_changed_S_to_R)) if ids_changed_S_to_R else ''
            })

    # After processing all variants, convert to DataFrame
    impact_df = pd.DataFrame(impact_results)

    # Order by Delta_Specificity ascending (most reduction first)
    impact_df_sorted = impact_df.sort_values(by='Delta_Specificity', ascending=True)

    # Save the impact metrics
    impact_file = os.path.join(output_dir, "variant_impact_metrics.tsv")
    impact_df_sorted.to_csv(impact_file, sep='\t', index=False)
    logging.info(f"Saved variant impact metrics to '{impact_file}'.")

    return impact_df_sorted  # Properly placed outside the loop and inside the function

def save_baseline_results(predictions_df, metrics, output_dir):
    """Saves baseline predictions and metrics."""
    baseline_dir = os.path.join(output_dir, "baseline")
    os.makedirs(baseline_dir, exist_ok=True)

    # Save predictions
    predictions_file = os.path.join(baseline_dir, "baseline_predictions.tsv")
    predictions_df.to_csv(predictions_file, sep='\t', index=False)
    logging.info(f"Saved baseline predictions to '{predictions_file}'.")

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_file = os.path.join(baseline_dir, "baseline_metrics.tsv")
    metrics_df.to_csv(metrics_file, sep='\t', index=False)
    logging.info(f"Saved baseline metrics to '{metrics_file}'.")

    # Save confusion matrix
    confusion_matrix = pd.DataFrame({
        'Actual': ['R', 'S'],
        'Predicted_R': [metrics['TP'], metrics['FP']],
        'Predicted_S': [metrics['FN'], metrics['TN']]
    })
    confusion_matrix_file = os.path.join(baseline_dir, "confusion_matrix.tsv")
    confusion_matrix.to_csv(confusion_matrix_file, sep='\t', index=False)
    logging.info(f"Saved baseline confusion matrix to '{confusion_matrix_file}'.")

def main():
    args = parse_arguments()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Setup logging with debug flag
    setup_logging(args.output_dir, debug=args.debug)
    logging.info("Starting Variant Impact Analysis.")

    # Load variant lists
    list1_df = load_list1_variants(args.list1_variants_file)
    list2_df = load_list2_variants(args.list2_variants_file)

    # Load ID list
    id_df = load_id_list(args.id_list_file)

    # Create combined variant mapping
    variant_map = create_variant_mapping(list1_df, list2_df)

    # Load all SNP files
    snp_dict, missing_snp = load_all_snp_files(id_df, args.snp_dir, args.additional_snp_dir)

    # Load all Indel files
    indel_dict, missing_indel = load_all_indel_files(id_df, args.indel_dir, args.additional_indel_dir)

    # Exclude IDs that do not have both SNP and Indel files
    included_ids = set(snp_dict.keys()).intersection(indel_dict.keys())
    excluded_ids = len(id_df) - len(included_ids)
    logging.info(f"Total IDs included (have both SNP and Indel files): {len(included_ids)}")
    logging.info(f"Total IDs excluded (do not have both SNP and Indel files): {excluded_ids}")

    # Filter snp_dict and indel_dict to include only IDs with both files
    snp_dict = {id: snp_set for id, snp_set in snp_dict.items() if id in included_ids}
    indel_dict = {id: indel_set for id, indel_set in indel_dict.items() if id in included_ids}

    if not included_ids:
        logging.error("No IDs have both SNP and Indel files. Exiting.")
        sys.exit(1)

    # Generate baseline predictions using only List1
    list1_set = set(zip(list1_df['position'], list1_df['ref'], list1_df['alt']))
    baseline_predictions_df, baseline_metrics = generate_baseline_predictions(id_df, snp_dict, indel_dict, list1_set)

    # Save baseline results
    save_baseline_results(baseline_predictions_df, baseline_metrics, args.output_dir)

    # Analyze impact of each List2 variant
    impact_df_sorted = analyze_variant_impact(list2_df, list1_set, snp_dict, indel_dict, baseline_metrics, id_df, args.output_dir)

    logging.info("Variant Impact Analysis completed successfully.")
    print(f"Total IDs processed: {len(included_ids)}")
    print(f"Total IDs excluded due to missing or unreadable SNP files: {missing_snp}")
    print(f"Total IDs excluded due to missing or unreadable Indel files: {missing_indel}")

if __name__ == "__main__":
    main()
