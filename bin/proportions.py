import numpy as np
from scipy import stats

def calculate_p_values(list1, list2):
  # Calculate the proportion of each category in each list
  proportions1 = [list1.count(i) / len(list1) for i in set(list1)]
  proportions2 = [list2.count(i) / len(list2) for i in set(list2)]

  # Calculate the number of permutations
  n_permutations = len(list1) * len(list2)

  # Create a combined list of proportions
  combined_proportions = proportions1 + proportions2

  # Calculate the observed difference in means
  observed_diff = np.mean(proportions1) - np.mean(proportions2)

  # Use a permutation test to calculate the p-values
  p_values = []
  for i in range(n_permutations):
    # Permute the combined proportions
    permuted_proportions = np.random.permutation(combined_proportions)

    # Calculate the difference in means for the permuted proportions
    permuted_diff = np.mean(permuted_proportions[:len(proportions1)]) - np.mean(permuted_proportions[len(proportions1):])

    # Calculate the p-value for the permuted proportions
    if permuted_diff >= observed_diff:
      p_values.append(1)
    else:
      p_values.append(0)

  # Calculate the average p-value
  p_value = np.mean(p_values)

  # Correct the p-value for multiple comparisons
  corrected_p_value = stats.multipletests(p_value, method='bonferroni')[1]

  # Calculate the false discovery rate
  fdr = stats.multipletests(p_value, method='fdr_bh')[1]

  return p_value, corrected_p_value, fdr

# Example:

list1 = ['a', 'b', 'a', 'c', 'c', 'd', 'e', 'e', 'e', 'f']
list2 = ['a', 'b', 'b', 'c', 'c', 'c', 'd', 'd', 'd', 'd']
p_value, corrected_p_value, fdr = calculate_p_values(list1, list2)

print("P-value:", p_value)
print("Corrected p-value:", corrected_p_value)
print("False discovery rate:", fdr)