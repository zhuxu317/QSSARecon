import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator

# Define line styles, colors, and markers
line_arr = ('-', '--', '-.', ':')
color_arr = ('k', 'b', 'm', 'y', 'g', 'c')
label_size = 10
font_size = 10
tick_size = 10
# Directory containing the CSV files
data_dir = "data/case_NH3_counterflow/ILDM"

# Create a figure with 2 subplots in a row
fig, axes = plt.subplots(1, 2, figsize=(5, 2.5))  # 2 subplots in 1 row

# Plot for OH, O, and H on the first subplot
species_to_manipulate_1 = ['OH', 'O', 'H']
species_color_map_1 = {species: color for species, color in zip(species_to_manipulate_1, color_arr)}
text_positions_1 = {
    'OH': (3e-4, 2.4e-3),
    'H': (7e-4, 0.8e-3),
    'O': (7e-4, 1.4e-4)
}

# Plot for NO and NO2 on the second subplot
species_to_manipulate_2 = ['NO', 'NO2']
species_color_map_2 = {species: color for species, color in zip(species_to_manipulate_2, color_arr)}
text_positions_2 = {
    'NO': (2e-2, 3e-5),
    'NO2': (9e-3, 3e-6)
}

fig_name = "figs/CEQ/ILDM_combined.pdf"

# Iterate over CSV files in the directory
for filename in os.listdir(data_dir):
    if filename.endswith(".csv"):
        file_path = os.path.join(data_dir, filename)
        species_name = filename.split('_')[0]

        # Read CSV data
        df = pd.read_csv(file_path)

        # Plot on the first subplot for OH, O, H
        if species_name in species_to_manipulate_1:
            ax1 = axes[0]  # First subplot (OH, O, H)
            ax1.plot(df['Time'], df[species_name], label=f'{species_name} - {filename}',
                     linestyle=line_arr[0], color=species_color_map_1.get(species_name, 'k'))
            # Add text annotation
            x_pos, y_pos = text_positions_1[species_name]
            ax1.text(x_pos, y_pos, f'{species_name}', color=species_color_map_1[species_name], fontsize=label_size)

        # Plot on the second subplot for NO, NO2
        elif species_name in species_to_manipulate_2:
            ax2 = axes[1]  # Second subplot (NO, NO2)
            ax2.plot(df['Time'], df[species_name], label=f'{species_name} - {filename}',
                     linestyle=line_arr[0], color=species_color_map_2.get(species_name, 'k'))
            # Add text annotation
            x_pos, y_pos = text_positions_2[species_name]
            ax2.text(x_pos, y_pos, species_name, color=species_color_map_2[species_name], fontsize=label_size)

# Add vertical line for residence time in both subplots
residence_time = 1/87

# Vertical line and text on the first subplot
axes[0].axvline(x=residence_time, color='red', linestyle='--')
axes[0].text(3e-4, 1e-5, '$t_{res}$', color='red', rotation=0, va='bottom', fontsize=label_size)

# Vertical line and text on the second subplot
axes[1].axvline(x=residence_time, color='red', linestyle='--')
axes[1].text(2.5e-4, 2e-9, '$t_{res}$', color='red', rotation=0, va='bottom', fontsize=label_size)

# Customize both subplots (logarithmic scale, labels, and tick formatting)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('t (s)', fontsize=font_size)
ax1.set_ylabel('X', fontsize=font_size)
ax1.set_xlim([1e-13, 1e-1])  # Set x-axis range for the first subplot
ax1.set_ylim([1e-6, 1e-1])   # Set y-axis range for the first subplot

# Reduce the number of ticks on both axes for ax1
ax1.xaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on x-axis
ax1.yaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on y-axis
ax1.xaxis.set_minor_locator(NullLocator())  # Disable minor ticks on x-axis
ax1.yaxis.set_minor_locator(NullLocator())  # Disable minor ticks on y-axis

# Increase font size of tick labels for ax1
ax1.tick_params(axis='both', which='major', labelsize=tick_size)

# Customize second subplot
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('t (s)', fontsize=font_size)
ax2.set_xlim([1e-13, 1e0])  # Set x-axis range for the second subplot
ax2.set_ylim([1e-9, 1e-2])   # Set y-axis range for the second subplot

# Reduce the number of ticks on both axes for ax2
ax2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on x-axis
ax2.yaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on y-axis
ax2.xaxis.set_minor_locator(NullLocator())  # Disable minor ticks on x-axis
ax2.yaxis.set_minor_locator(NullLocator())  # Disable minor ticks on y-axis

# Increase font size of tick labels for ax2
ax2.tick_params(axis='both', which='major', labelsize=tick_size)

# Adjust layout and reduce the spacing between the subplots
plt.tight_layout()
# plt.subplots_adjust(wspace=0.2)  # Adjust wspace to control horizontal gap

# Show the plot
plt.show()
plt.ioff()

# Save the figure with tight bounding box to remove edges
plt.savefig(fig_name, dpi=300, bbox_inches='tight')