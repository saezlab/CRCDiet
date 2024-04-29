"""import numpy as np
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
print("False discovery rate:", fdr)"""

# given two lists where items are repeated, calculate the p-value of the difference in proportions of each item in the two lists in python, write the code below
#   - calculate the proportion of each item in each list
#  - calculate the number of permutations
# - create a combined list of proportions
#   
import numpy as np
# Create frequency tables for each list
def create_proportion_table(data):
    categories = set(data)
    freq_table = {}
    for category in categories:
        freq_table[category] = data.count(category)/len(data)
        print(data.count(category), len(data))
    return freq_table
"""# Define your two lists with categories
list1 = ["Category1", "Category1", "Category1", "Category2", "Category1", "Category3", "Category2", "Category1", "Category4"]
list2 = ["Category1", "Category2", "Category2", "Category4", "Category4","Category4","Category4", "Category3", "Category3", "Category4", "Category4"]



freq_table1 = create_proportion_table(list1)
freq_table2 = create_proportion_table(list2)


# Create a list of all unique categories
unique_categories = set(freq_table1.keys()).union(set(freq_table2.keys()))

# Initialize a dictionary to store p-values for each category
p_values = {}

# Define the number of permutations
num_permutations = 100

# Perform the permutation test for each category individually
for category in unique_categories:
    print(category, freq_table1[category], freq_table2[category])
    observed_diff = freq_table1.get(category, 0) - freq_table2.get(category, 0)
    print(observed_diff)

    all_data = list1 + list2
    num_obs = len(list1)
    num_perm = 0
    num_of_rare_or_rarer = 0.0
    for _ in range(num_permutations):
        np.random.shuffle(all_data)
        # print(all_data)
        perm_list1 = all_data[:num_obs]
        perm_list2 = all_data[num_obs:]

        perm_diff = perm_list1.count(category)/len(perm_list1) - perm_list2.count(category)/len(perm_list2)
        if observed_diff <0.0 and perm_diff<observed_diff:

            num_of_rare_or_rarer += 1
            print("First ", category,  observed_diff, perm_diff)
        elif observed_diff >=0.0 and perm_diff>observed_diff:
            num_of_rare_or_rarer +=1
            print("Second ", category,  observed_diff, perm_diff)

    # print(num_of_rare_or_rarer)
    p_value = (num_of_rare_or_rarer + 1) / (num_permutations + 1)
    p_values[category] = p_value

# Set a significance level (e.g., 0.05)
significance_level = 0.05

# Report whether the change in proportion is statistically significant for each category
for category, p in p_values.items():
    if p < significance_level:
        print(f"Category: {category} - p-value: {p} (Significant change)")
    else:
        print(f"Category: {category} - p-value: {p} (No significant change)")"""

def cal_significance_proportions(list1, list2, num_permutations=10000):
    freq_table1 = create_proportion_table(list1)
    freq_table2 = create_proportion_table(list2)
    # Create a list of all unique categories
    unique_categories = set(freq_table1.keys()).union(set(freq_table2.keys()))

    # Initialize a dictionary to store p-values for each category
    p_values = {}


    # Perform the permutation test for each category individually
    for category in unique_categories:
        # print(category, freq_table1[category], freq_table2[category])
        observed_diff = freq_table1.get(category, 0) - freq_table2.get(category, 0)
        # print(observed_diff)

        all_data = list1 + list2
        num_obs = len(list1)
        num_perm = 0
        num_of_rare_or_rarer = 0.0
        for _ in range(num_permutations):
            np.random.shuffle(all_data)
            # print(all_data)
            perm_list1 = all_data[:num_obs]
            perm_list2 = all_data[num_obs:]

            perm_diff = perm_list1.count(category)/len(perm_list1) - perm_list2.count(category)/len(perm_list2)
            """print("First proportion, second p", perm_list1.count(category)/len(perm_list1), perm_list2.count(category)/len(perm_list2))
            print("Cat, perm_diff, obs diff", category, perm_diff, observed_diff)
            print("observed_diff", observed_diff)
            print("Perm. diff", perm_diff)"""
            if observed_diff <0.0 and perm_diff<observed_diff:

                num_of_rare_or_rarer += 1
                # print("First ", category,  observed_diff, perm_diff)
            elif observed_diff >=0.0 and perm_diff>observed_diff:
                num_of_rare_or_rarer +=1
                # print("Second ", category,  observed_diff, perm_diff)

        # print(num_of_rare_or_rarer)
        p_value = (num_of_rare_or_rarer + 1) / (num_permutations + 1)
        p_values[category] = p_value

    # Set a significance level (e.g., 0.05)
    significance_level = 0.05

    # Report whether the change in proportion is statistically significant for each category
    for category, p in p_values.items():
        if p < significance_level:
            print(f"Category: {category} - p-value: {p} (Significant change)")
        else:
            print(f"Category: {category} - p-value: {p} (No significant change)")
            

list1 = ["Category1", "Category1", "Category1", "Category2", "Category1", "Category3", "Category2", "Category1", "Category4"]
list2 = ["Category1", "Category2", "Category2", "Category4", "Category4","Category4","Category4", "Category3", "Category3", "Category4", "Category4"]

cal_significance_proportions(list1, list2, num_permutations=10000)
