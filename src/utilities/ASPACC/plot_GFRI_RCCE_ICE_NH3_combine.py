import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator

method_path = "figs/GFRI/data/"
# Define target cases
target_case_assemble = [
    'N_CF_1',
    'N_CF_2',
    'N_CF_3',
    'N_CF_4',
    'N_CF_5',
    'N_CF_6',
    'N_CF_7',
    'N_CF_8',
]
target_name_assemble = [
    'N_CF_1',
    'N_CF_2',
    'N_CF_3',
    'N_CF_4',
    'N_CF_5',
    'N_CF_6',
    'N_CF_7',
    'N_CF_8',
]

x_range_assemble = [
    [18, 27],
    [18, 35],
    [17, 26],
    [17, 26],
    [18, 35],
    [22, 31],
    [22, 31],
    [22, 31],
]

case_name = "case_NH3_counterflow"

line_arr = ('-','--','-.',':')
color_arr = ('k','r','b','y','g', 'c','m')
# color_arr = ('k','r','b','g', 'm')

symbol_arr = ('s','o','v','^','*')
# #############################FIGURE OVERALL OVER T#########################
# fig, axs = plt.subplots(3, 3, figsize=(8, 7), sharex=True)

# data_cols = ['OH', 'O', 'H', 'N', 'NO', 'NO2', 'CEM', 'MF', 'Qdot']

# for j, target_case in enumerate(target_case_assemble):
#     # Define paths
#     original_data_path = os.path.join("data", case_name, target_case) + ".csv"
#     rcce_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")
#     gfri_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")
#     ice_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")

#     # Read data
#     df_original = pd.read_csv(original_data_path)
#     df_rcce = pd.read_csv(rcce_data_path)
#     df_gfri = pd.read_csv(gfri_data_path)
#     df_ice = pd.read_csv(ice_data_path)

#     T = df_original['T']  # Temperature data for x-axis

#     # Plot each variable in its corresponding subplot
#     for i, col in enumerate(data_cols):
#         row, col_num = divmod(i, 3)
#         ax = axs[row, col_num]
#         # Scatter plot for original data
#         ax.scatter(T[::3], df_original[col][::3], label=f'{target_case}', s=10, alpha=0.8, marker=symbol_arr[j % len(symbol_arr)], color=color_arr[j % len(color_arr)])
#         # Line plot for RCCE, GFRI, and ICE for each case
#         # ax.plot(T, df_rcce[col], label=f'{target_case} - RCCE', linestyle='--', linewidth=1, color=colors[i])
#         ax.plot(T, df_gfri[col], linestyle=line_arr[j % len(line_arr)], linewidth=0.8, color=color_arr[j % len(color_arr)])
#         # ax.plot(T, df_ice[col], label=f'{target_case} - ICE', linestyle='-', linewidth=1, color=colors[i])
#         # Add labels
#         ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#         ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#         # ax.xaxis.get_major_formatter().set_powerlimits((-2, 2))
#         ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
#         ax.set_ylabel(col, fontsize=10)
#         ax.set_xlim([1000, 2500])
#         if row == 2:
#             ax.set_xlabel('T(K)', fontsize=10)

# for j in range(len(data_cols), 9):
#     fig.delaxes(axs[j // 3, j % 3])
# plt.tight_layout()

# # Save the figure
# fig_dir ="figs/NH3"
# os.makedirs(fig_dir, exist_ok=True)
# plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_Temperature.png"), dpi=300)


# ###################### Plot FIGURE RECONSTRUCT SPACE ######################
fig, axs = plt.subplots(3, len(target_case_assemble), figsize=(3 * len(target_case_assemble), 8))

# Ensure axs is always 2D even if there's only one column
if len(target_case_assemble) == 1:
    axs = axs[:, np.newaxis]

cols_1 = ['T', 'NH3', 'H2', 'O2', 'N2']
cols_2 = ['MF', 'OH', 'O','H', 'NO']
cols_3 = ['CEM', 'Qdot']

for j, target_case in enumerate(target_case_assemble):
    original_data_path = os.path.join("data", case_name, target_case + ".csv")
    gfri_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")

    # Read data with error handling for file paths
    try:
        df_original = pd.read_csv(original_data_path)
        df_gfri = pd.read_csv(gfri_data_path)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        continue

    grid = df_original['grid'] * 1000
    # Plot cols_1 (first row)
    ax1 = axs[0, j]
    ax2 = ax1.twinx()
    for i, col in enumerate(cols_1):
        if col == 'T':
            ax2.plot(
                grid, 
                df_original[col], 
                linestyle=line_arr[i % len(line_arr)], 
                linewidth=2, 
                color='red', 
                label='T'
            )
        else:
            ax1.plot(
                grid, 
                df_gfri[col], 
                linestyle=line_arr[(i+1) % len(line_arr)], 
                linewidth=2, 
                color=color_arr[(i+1) % len(color_arr)], 
                label=f'{col}'
            )
            # Optional: Add scatter if needed
            # ax1.scatter(
            #     grid, 
            #     df_original[col], 
            #     linestyle=line_arr[(i+1) % len(line_arr)], 
            #     linewidth=2, 
            #     color=color_arr[(i+1) % len(color_arr)], 
            #     marker=symbol_arr[j % len(symbol_arr)]
            # )

    # Set axis limits
    ax1.set_xlim(x_range_assemble[j])
    ax2.set_xlim(x_range_assemble[j])
    ax2.set_ylim([0, 2500])
    ax1.set_ylim([0, 0.9])

    # Configure y-labels and ticks
    if j == 0:
        ax1.set_ylabel('X', fontsize=10, color='black')
        ax1.tick_params(axis='y', labelleft=True, labelright=False, direction='in')
        ax2.tick_params(axis='y', labelleft=False, labelright=False)
    elif j == len(target_case_assemble) - 1:
        ax2.set_ylabel('T(K)', fontsize=10, color='red')
        
        # Collect handles and labels from both ax2 and ax1
        lines_ax2, labels_ax2 = ax2.get_legend_handles_labels()
        lines_ax1, labels_ax1 = ax1.get_legend_handles_labels()
        
        # Combine them with ax2's labels first
        combined_lines = lines_ax2 + lines_ax1
        combined_labels = labels_ax2 + labels_ax1
        
        # Create the combined legend
        legend = ax1.legend(
            combined_lines, 
            combined_labels, 
            loc='center left', 
            bbox_to_anchor=(0.55, 0.55), 
            frameon=False, 
            fontsize=12
        )
        
        # Set legend text colors to match the line colors
        for text, line in zip(legend.get_texts(), legend.get_lines()):
            text.set_color(line.get_color())
        
        # Configure tick parameters
        ax2.tick_params(axis='y', labelleft=False, colors='red', direction='in', labelright=True)
        ax1.tick_params(axis='y', labelleft=False, labelright=False)
    else:
        ax2.tick_params(axis='y', labelleft=False, labelright=False)
        ax1.tick_params(axis='y', labelleft=False, labelright=False)

    # Configure x-axis ticks and titles
    ax1.tick_params(axis='x', labelbottom=False)
    ax1.set_title(f'{target_name_assemble[j]}', fontsize=10)
    
    # ###################### Plot cols_2 (second row) ######################
    ax1_2 = axs[1, j]
    ax2_2 = ax1_2.twinx()
    for i, col in enumerate(cols_2):
        marker = symbol_arr[i % len(symbol_arr)]
        if col == 'MF':
            ax1_2.scatter(
                grid[::3], 
                df_original[col][::3], 
                s=10, 
                alpha=0.8, 
                marker=marker, 
                color=color_arr[i % len(color_arr)]
            ) 
            ax1_2.plot(
                grid, 
                df_gfri[col], 
                linestyle=line_arr[i % len(line_arr)], 
                linewidth=1, 
                color='black', 
                label='MF'
            )
        else:
            ax2_2.scatter(
                grid[::3], 
                df_original[col][::3], 
                s=10, 
                alpha=0.8, 
                marker=marker, 
                color=color_arr[i % len(color_arr)]
            )
            ax2_2.plot(
                grid, 
                df_gfri[col], 
                linestyle=line_arr[i % len(line_arr)], 
                linewidth=1.2, 
                color=color_arr[i % len(color_arr)], 
                label=f'{col}'
            )
    
    # Set axis limits and scales
    ax1_2.set_xlim(x_range_assemble[j])
    ax2_2.set_xlim(x_range_assemble[j])
    ax1_2.set_ylim([0, 1])
    ax2_2.set_ylim([1e-10, 1e-2])
    ax2_2.set_yscale('log')  # Set y-axis to log scale for ax2
    
    # Configure y-labels and ticks for cols_2
    if j == 0:
        ax1_2.set_ylabel('MF', fontsize=10, color='black')
        ax1_2.tick_params(axis='y', labelleft=True, colors='black', direction='in', labelright=False)
        ax2_2.tick_params(axis='y', labelleft=False, labelright=False, direction='in')
    elif j == len(target_case_assemble) - 1:
        ax2_2.set_ylabel('X', fontsize=10)
        
        # Collect handles and labels from both ax2_2 and ax1_2
        lines_ax2_2, labels_ax2_2 = ax2_2.get_legend_handles_labels()
        lines_ax1_2, labels_ax1_2 = ax1_2.get_legend_handles_labels()
        
        # Combine them with ax2_2's labels first
        combined_lines_2 = lines_ax2_2 + lines_ax1_2
        combined_labels_2 = labels_ax2_2 + labels_ax1_2
        
        # Create the combined legend for cols_2
        legend_2 = ax1_2.legend(
            combined_lines_2, 
            combined_labels_2, 
            loc='center left', 
            bbox_to_anchor=(0.55, 0.55), 
            frameon=False, 
            fontsize=12
        )
        
        # Set legend text colors to match the line colors
        for text, line in zip(legend_2.get_texts(), legend_2.get_lines()):
            text.set_color(line.get_color())
        
        # Configure tick parameters
        ax2_2.tick_params(axis='y', labelleft=False, direction='in', labelright=True)
        ax1_2.tick_params(axis='y', labelleft=False, labelright=False)
    else:
        ax1_2.tick_params(axis='y', labelleft=False, labelright=False)
        ax2_2.tick_params(axis='y', labelleft=False, labelright=False)
    
    # Configure x-axis ticks
    ax1_2.tick_params(axis='x', labelbottom=False)
    # Plot cols_3 (third row)
    ax1, ax2 = axs[2, j], axs[2, j].twinx()

    # Lists to store the legend handles and labels
    legend_handles = []
    legend_labels = []

    for i, col in enumerate(cols_3):
        marker = symbol_arr[i % len(symbol_arr)]
        
        # Convert the pandas Series to numpy array before indexing
        grid_values = grid  # Assuming grid is already a numpy array, if not, convert it to one
        df_original_col_values = df_original[col].to_numpy()  # Convert pandas series to numpy array
        df_gfri_col_values = df_gfri[col].to_numpy()  # Convert pandas series to numpy array

        if col == 'CEM':
            ax1.plot(grid_values, df_gfri_col_values, linestyle=line_arr[i % len(line_arr)], 
                    linewidth=1, color='black', label='CEM')
            ax1.scatter(grid_values[::2], df_original_col_values[::2], s=10, alpha=0.8, 
                        marker=marker, color=color_arr[i % len(color_arr)]) 
            # Add handle for the legend
            legend_handles.append(ax1.lines[-1])  # Last plot line added to ax1
            legend_labels.append('CEM')
        else:
            ax2.scatter(grid_values[::2], df_gfri_col_values[::2], s=10, alpha=0.8, 
                        marker=marker, color=color_arr[i % len(color_arr)])
            ax2.plot(grid_values, df_original_col_values, linestyle=line_arr[i % len(line_arr)], 
                    linewidth=1.2, color=color_arr[i % len(color_arr)], label=f'HRR')
            # Add handle for the legend
            legend_handles.append(ax2.lines[-1])  # Last plot line added to ax2
            legend_labels.append(f'HRR')

    # Set axis limits
    ax1.set_xlim(x_range_assemble[j])
    ax2.set_xlim(x_range_assemble[j])
    ax1.set_ylim([-5, 2])
    ax2.set_ylim([0, 6e8])

    # Axis labeling and tick parameters
    if j == 0:
        ax1.set_ylabel('CEM', fontsize=10, color='black')
        ax1.tick_params(axis='y', labelleft=True, colors='black', direction='in', labelright=False)
        ax2.tick_params(axis='y', labelleft=False, labelright=False, direction='in')
    elif j == len(target_case_assemble) - 1:
        ax2.set_ylabel('HRR(J/kg/s)', fontsize=10, color='red')

        # Combine the legend for both axes
        legend_2 = ax2.legend(handles=legend_handles, labels=legend_labels, loc='center left', bbox_to_anchor=(0.55, 0.55), frameon=False, fontsize=12)

        # Set legend text colors to match the line colors
        for text, line in zip(legend_2.get_texts(), legend_2.get_lines()):
            text.set_color(line.get_color())
        
        ax1.tick_params(axis='y', labelleft=False, labelright=False)
        ax2.tick_params(axis='y', labelleft=False, direction='in', labelright=True, colors='red')
    else:
        ax1.tick_params(axis='y', labelleft=False, labelright=False)
        ax2.tick_params(axis='y', labelleft=False, labelright=False)

    # Add x-labels for the third row
    ax1.set_xlabel('x (mm)', fontsize=10)
    ax1.tick_params(axis='x', labelbottom=True)



plt.tight_layout()
plt.show()
# Save the figure
fig_dir = "figs/NH3"
os.makedirs(fig_dir, exist_ok=True)
plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_expanded.pdf"), dpi=300, bbox_inches='tight')

#############################FIGURE OVERALL OVER T#########################
# fig, axs = plt.subplots(2, 2, figsize=(5, 4), sharex=True)
# data_cols = ['OH', 'CEM', 'MF', 'Qdot']

# for j, target_case in enumerate(target_case_assemble):
#     # Define paths
#     original_data_path = os.path.join("data", case_name, target_case) + ".csv"
#     rcce_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")
#     gfri_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")
#     ice_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")

#     # Read data
#     df_original = pd.read_csv(original_data_path)
#     df_rcce = pd.read_csv(rcce_data_path)
#     df_gfri = pd.read_csv(gfri_data_path)
#     df_ice = pd.read_csv(ice_data_path)

#     T = df_original['T']  # Temperature data for x-axis

#     # Plot each variable in its corresponding subplot
#     for i, col in enumerate(data_cols):
#         row, col_num = divmod(i, 2)
#         ax = axs[row, col_num]
#         # Scatter plot for original data
#         ax.scatter(T[::3], df_original[col][::3], label=f'{target_case}', s=10, alpha=0.8, color=color_arr[j % len(color_arr)], marker=symbol_arr[j % len(symbol_arr)])
#         ax.plot(T, df_gfri[col], linestyle=line_arr[j % len(line_arr)], linewidth=0.8, color=color_arr[j % len(color_arr)])
#         if col == "Qdot":
#             ax.set_ylabel("HRR(J/kg/s)", fontsize=10)
#         else:
#             ax.set_ylabel(col, fontsize=10)
#         ax.set_xlim([1000, 2500])
#         ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#         ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#         ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
#         if row == 1:
#             ax.set_xlabel('T(K)', fontsize=10)

#     # Add target_case text outside the figure
    
#     #     text_positions = [
#     #     (0.35, 0.80),  # Position for the 1st text
#     #     (0.365, 0.84),  # Position for the 2nd text
#     #     (0.38, 0.88),  # Position for the 3rd text
#     # ]
#     text_positions = [
#         (0.36, 0.88),  # Position for the 1st text
#         (0.34, 0.83),  # Position for the 2nd text
#         (0.325, 0.78),  # Position for the 3rd text
#         (0.28, 0.70),  # Position for the 4th text
#         (0.25, 0.66),  # Position for the 5th text
#     ]

#     # Ensure we are within bounds of the positions array
#     if j < len(text_positions):
#         fig.text(
#             text_positions[j][0],  # X position from predefined positions
#             text_positions[j][1],  # Y position from predefined positions
#             target_name_assemble[j],  # Text to display
#             ha='center', va='center', fontsize=8, 
#             color=color_arr[j % len(color_arr)], 
#             transform=fig.transFigure
#         )

# # Remove any unused subplots if less than 4 data columns
# for j in range(len(data_cols), 4):
#     fig.delaxes(axs[j // 2, j % 2])

# plt.tight_layout()
# fig_dir = "figs/NH3"
# os.makedirs(fig_dir, exist_ok=True)
# plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_Result2.pdf"), dpi=300)

# ###################### Plot FIGURE DECREASING MOLE FRACTION and ERROR#####################

# # Initialize a Series to accumulate species data
# all_species_data = pd.Series(dtype=float)

# # Create a figure and primary axis
# fig, ax1 = plt.subplots(figsize=(6, 4))  # Increased figsize for clarity

# # Loop to accumulate maximum species data across all cases
# for j, target_case in enumerate(target_case_assemble):
#     original_data_path = os.path.join("data", case_name, target_case) + ".csv"
#     gfri_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")

#     # Read data with error handling
#     try:
#         df_original = pd.read_csv(original_data_path)
#         df_gfri = pd.read_csv(gfri_data_path)
#     except FileNotFoundError as e:
#         print(f"File not found: {e}")
#         continue

#     # Keep only the specified columns
#     keep_cols = ['Qdot','MF', 'OH', 'NO', 'H', 'O', 'NH2', 'N2O', 'NO2', 'NH', 'N2H2', 'N']
#     species_data = df_original[keep_cols]
#     gfri_species_data = df_gfri[keep_cols]

#     # Accumulate maximum species data
#     all_species_data = all_species_data.add(species_data.max(), fill_value=0)

# # Calculate the average mole fraction for each species across all cases
# average_mole_fractions = all_species_data / len(target_case_assemble)

# # Sort species by decreasing average mole fraction
# sorted_species = average_mole_fractions.sort_values(ascending=False)

# # Create a separate copy for the bar plot and set 'Qdot' and 'MF' to 0
# sorted_species_bar = sorted_species.copy()
# sorted_species_bar.loc['Qdot'] = 0
# sorted_species_bar.loc['MF'] = 0

# # Rename 'Qdot' to 'HRR' and 'MF' to 'X' for labeling in bar plot
# sorted_species_bar.rename(index={'Qdot': 'HRR'}, inplace=True)

# # Bar plot for the average mole fraction using sorted_species_bar
# ax1.bar(
#     range(len(sorted_species_bar)),  # Numeric x positions
#     sorted_species_bar.values,
#     edgecolor='black',
#     facecolor='none',
#     hatch='//',
#     linewidth=0.5,
#     label='Average Mole Fraction'
# )
# ax1.set_ylabel("Average Mole Fraction", fontsize=12)
# ax1.set_yscale("log")
# ax1.set_ylim(1e-6, 1e-2)  # Set y-axis range
# ax1.tick_params(axis='x', rotation=0, labelsize=10)
# ax1.tick_params(axis='y', labelsize=10)

# # Set x-axis labels: 'HRR' and 'X' along with others if present
# new_labels = sorted_species_bar.index.tolist()
# ax1.set_xticks(range(len(sorted_species_bar)))
# ax1.set_xticklabels(new_labels, rotation=0, fontsize=10)

# # Create a secondary y-axis once, outside the loop
# ax2 = ax1.twinx()

# # Loop to plot the max ratio error for each target case
# for j, target_case in enumerate(target_case_assemble):
#     original_data_path = os.path.join("data", case_name, target_case) + ".csv"
#     gfri_data_path = os.path.join(method_path, case_name, target_case, "predicted_X.csv")

#     # Read data with error handling
#     try:
#         df_original = pd.read_csv(original_data_path)
#         df_gfri = pd.read_csv(gfri_data_path)
#     except FileNotFoundError as e:
#         print(f"File not found: {e}")
#         continue

#     species_data = df_original[keep_cols]
#     gfri_species_data = df_gfri[keep_cols]

#     # Calculate the max ratio error
#     # Handle potential division by zero by replacing infinities and NaNs
#     with np.errstate(divide='ignore', invalid='ignore'):
#         max_ratio_error = ((gfri_species_data.max() - species_data.max()) / species_data.max()) * 100
#         max_ratio_error.replace([np.inf, -np.inf], np.nan, inplace=True)
#         max_ratio_error.fillna(0, inplace=True)
    
#     # Map 'HRR' and 'X' back to 'Qdot' and 'MF' for the line plots
#     max_ratio_error_mapped = max_ratio_error.copy()
#     max_ratio_error_mapped = max_ratio_error_mapped.rename(index={'HRR': 'Qdot', 'X': 'MF'})

#     # Plot on ax2 using numeric x positions
#     ax2.plot(
#         range(len(sorted_species_bar)),  # Numeric x positions
#         max_ratio_error_mapped.values,
#         marker=symbol_arr[j % len(symbol_arr)],
#         linestyle=line_arr[j % len(line_arr)],
#         color=color_arr[j % len(color_arr)],  # Use color for the edge of the marker
#         markerfacecolor='none',  # Hollow marker
#         label=f'{target_name_assemble[j]}',
#         linewidth=1.0,
#         markersize=6
#     )

# # Set labels and limits for ax2
# ax2.set_ylabel("Peak Value Reconstruction Error(%)", fontsize=12)
# ax2.set_ylim(-150, 500)  # Adjust the range as needed
# ax2.tick_params(axis='y', labelsize=10)

# # Add ±20% and ±50% error lines with annotations
# error_lines = [
#     (20, 'red', '20%'),
#     (-20, 'red', ''),
#     (50, 'blue', '50%'),
#     (-50, 'blue', ''),
# ]

# for y, color, label_text in error_lines:
#     ax2.axhline(y=y, color=color, linestyle='--', linewidth=1)
#     ax2.text(len(sorted_species_bar) - 1.5, y, label_text, color=color, fontsize=10, 
#              verticalalignment='bottom' if y > 0 else 'top')

# # Set the tick label fontweight and size for ax2 to match ax1
# for label in ax2.get_yticklabels():
#     label.set_fontweight('normal')
#     label.set_fontsize(10)  # Ensure the font size matches ax1

# # Ensure ax1 tick labels also have normal font weight and size
# for label in ax1.get_yticklabels():
#     label.set_fontweight('normal')
#     label.set_fontsize(10)

# # Add a legend for ax2
# legend = ax2.legend(loc='upper right', frameon=False, fontsize=10)
# # Set legend text colors to match the line colors
# for text, line in zip(legend.get_texts(), legend.get_lines()):
#     text.set_color(line.get_color())

# # Optionally, add a legend for the bar plot (ax1)
# bar_legend = ax1.legend(loc='upper left', frameon=False, fontsize=10)
# bar_legend.get_texts()[0].set_color('black')  # Assuming only one label

# # Adjust layout to prevent overlap
# plt.tight_layout()

# # Save the figure
# fig_dir = "figs/NH3"
# os.makedirs(fig_dir, exist_ok=True)
# plt.savefig(os.path.join(fig_dir, "Average_Mole_Fraction_with_Max_Ratio_Error.pdf"), dpi=300)
# plt.show()
# plt.close(fig)
