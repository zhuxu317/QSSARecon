import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Load the input and output CSV files
input_file = 'data/1D/test/main_species_all.csv'  # Replace with your actual path
output_file = 'data/1D/test/minor_species_all.csv'  # Replace with your actual path

# Read input and output data
input_data = pd.read_csv(input_file)
output_data = pd.read_csv(output_file)

# Get input and output species from the headers of the CSV
input_species = input_data.columns.tolist()
output_species = output_data.columns.tolist()

# Step 2: Select the original input and output from the 20,000th row
original_input = input_data.iloc[20000].values  # Original input at row 20000
original_output = output_data.iloc[20000].values  # Original output at row 20000

# Step 3: Filter data where all input values are within ±5% of the original input
perturbation_rate = 0.05
lower_bound = original_input * (1 - perturbation_rate)
upper_bound = original_input * (1 + perturbation_rate)

# Filter the input data where each column value is within 5% of the original input
mask = np.all((input_data >= lower_bound) & (input_data <= upper_bound), axis=1)
filtered_input_data = input_data[mask]
filtered_output_data = output_data[mask]

# Step 4: Extract output data for specific species (OH, CH2O, NH, NH2) for filtered rows
output_selection = ["OH", "CH2O", "NH", "NH2"]
output_indices = [output_species.index(species) for species in output_selection]
selected_output_data = filtered_output_data.iloc[:, output_indices]

# Step 5: Plot violin plots for both selected input and output data

# Plot violin plot for filtered input data (all input species)
plt.figure(figsize=(10, 6))
sns.violinplot(data=filtered_input_data, scale='width')
plt.title("Violin Plot of Filtered Input Data within 5% of Original")
plt.xlabel("Input Species")
plt.ylabel("Value")
plt.xticks(ticks=np.arange(len(input_species)), labels=input_species, rotation=45)
plt.show()

# Plot violin plot for selected output data (OH, CH2O, NH, NH2)
plt.figure(figsize=(10, 6))
sns.violinplot(data=selected_output_data, scale='width')
plt.title("Violin Plot of Corresponding Output Data (OH, CH2O, NH, NH2)")
plt.xlabel("Output Species")
plt.ylabel("Concentration")
plt.xticks(ticks=np.arange(len(output_selection)), labels=output_selection)
plt.yscale('log')  # Use log scale for the y-axis
plt.show()

# Optionally, save the plots
plt.savefig('figs/sensitivity_analysis_output.png')

# Optionally, print the original input and filtered data
print("Original Input (Row 20000):", original_input)
print("Filtered Input Data within ±5% range:", filtered_input_data)
