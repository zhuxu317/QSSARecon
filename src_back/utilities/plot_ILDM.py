import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullLocator

# Define line styles, colors, and markers
line_arr = ('-', '--', '-.', ':')
color_arr = ('k','b','m','y','g', 'c')
# num_colors = 5
# viridis = plt.cm.get_cmap('jet', num_colors)
# color_arr = [viridis(i) for i in range(num_colors)]
label_size = 10
font_size = 10
tick_size = 10
# Directory containing the CSV files
data_dir = "data/case_NH3_counterflow/ILDM"

# Create a single figure
fig, ax = plt.subplots(figsize=(3.3, 2.8))  # Single plot

# Species to plot
species_to_manipulate = ['OH', 'O', 'H', 'NO', 'NO2']
species_color_map = {species: color for species, color in zip(species_to_manipulate, color_arr)}
text_positions = {
    'OH': (3e-4, 2.4e-3),
    'H': (7e-4, 0.8e-3),
    'O': (7e-4, 1.4e-4),
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

        # Plot species data if in the list
        if species_name in species_to_manipulate:
            ax.plot(df['Time'], df[species_name], label=f'{species_name} - {filename}',
                    linestyle=line_arr[0], color=species_color_map.get(species_name, 'k'))
            # Add text annotation
            x_pos, y_pos = text_positions[species_name]
            ax.text(x_pos, y_pos, species_name, color=species_color_map[species_name], fontsize=label_size)

# Add vertical line for residence time
residence_time = 1/87
ax.axvline(x=residence_time, color='red', linestyle='--')
ax.text(3e-4, 1e-7, '$t_{res}$', color='red', rotation=0, va='bottom', fontsize=label_size)

# Customize the plot (logarithmic scale, labels, and tick formatting)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('t (s)', fontsize=font_size)
ax.set_ylabel('X', fontsize=font_size)
ax.set_xlim([1e-13, 1e0])  # Set x-axis range
ax.set_ylim([1e-8, 1e-1])   # Set y-axis range

# Reduce the number of ticks on both axes
ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on x-axis
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=3))  # Reduce major ticks on y-axis
ax.xaxis.set_minor_locator(NullLocator())  # Disable minor ticks on x-axis
ax.yaxis.set_minor_locator(NullLocator())  # Disable minor ticks on y-axis

# Increase font size of tick labels
ax.tick_params(axis='both', which='major', labelsize=tick_size)

# Adjust layout
plt.tight_layout()

# Show the plot
plt.show()
plt.ioff()

# Save the figure with tight bounding box to remove edges
plt.savefig(fig_name, dpi=300, bbox_inches='tight')