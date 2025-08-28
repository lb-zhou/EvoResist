#!/usr/bin/env python3

## This script is used to:
## Perform a permutation test to evaluate whether the difference between the observed event count and 100 simulated event counts is statistically significant.

import numpy as np
from scipy.stats import anderson_ksamp
import pandas as pd

# Load observed data efficiently
# observed_summary_v3_filter.csv was generated using Table S5.
df = pd.read_csv("/work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/observed_summary_v3_filter.csv") 
df.columns = ["Event_number", "n", "Frequency"]

# Simulated dataset (frequency-based, no expansion)
simulated_data = {
    1: 196479048,
    2: 17414950,
    3: 1029117,
    4: 45963,
    5: 1665,
    6: 49,
    7:2
}

# Convert observed dataset into a dictionary
observed_data = dict(zip(df["Event_number"], df["n"]))

# Compute mean homoplastic mutations per sample (normalized for different sample sizes)
total_simulated_mutations = sum(k * v for k, v in simulated_data.items())
total_observed_mutations = sum(k * v for k, v in observed_data.items())

total_simulated_samples = sum(simulated_data.values())
total_observed_samples = sum(observed_data.values())

mean_simulated = total_simulated_mutations / total_simulated_samples
mean_observed = total_observed_mutations / total_observed_samples

# Convert dictionary to arrays (without expanding)
simulated_values = np.array(list(simulated_data.keys()))
simulated_weights = np.array(list(simulated_data.values()))

observed_values = np.array(list(observed_data.keys()))
observed_weights = np.array(list(observed_data.values()))

# Combine datasets for permutation test
combined_samples = np.concatenate([simulated_values, observed_values])

# Perform Anderson-Darling test (default method due to SciPy version)
observed_result = anderson_ksamp([simulated_values, observed_values])
observed_statistic = observed_result.statistic

# Permutation test (optimized)
num_permutations = 10000
permuted_stats = []

for i in range(num_permutations):
    # Randomly sample with replacement from the combined dataset
    perm_sample1 = np.random.choice(combined_samples, size=len(simulated_values), replace=True)
    perm_sample2 = np.random.choice(combined_samples, size=len(observed_values), replace=True)
    
    stat = anderson_ksamp([perm_sample1, perm_sample2]).statistic
    permuted_stats.append(stat)
    
    # Print progress every 100 iterations
    if (i + 1) % 100 == 0:
        print(f"Permutation {i + 1}/{num_permutations} completed...")

print("Permutation test finished!")

# Compute exact empirical p-value
p_value_exact = np.mean(np.array(permuted_stats) >= observed_statistic)

# Print results
print(f"Mean Homoplastic Mutations (Simulated): {mean_simulated:.4f}")
print(f"Mean Homoplastic Mutations (Observed): {mean_observed:.4f}")
print(f"Exact p-value (Monte Carlo Permutation Test): {p_value_exact:.6f}")

# Final sentence:
print(f"This simulation returned a skewed distribution that the number of homoplastic mutations "
      f"is significantly smaller than the observed ({mean_simulated:.4f} vs {mean_observed:.4f}, "
      f"P={p_value_exact:.6f}, Fig 2B).")
