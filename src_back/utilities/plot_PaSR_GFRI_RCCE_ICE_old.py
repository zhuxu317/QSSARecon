import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

inch = 2.54
x_range = [0, 0.2]
target_case = "DelftIII_MC_CEM"
original_data_path = os.path.join("data/PaSR", target_case) + ".csv"
rcce_data_path = os.path.join("figs/RCCE/data/PaSR", target_case, "predicted_Y.csv")
gfri_data_path = os.path.join("figs/GFRI/data/PaSR", target_case, "predicted_Y.csv")
# ice_data_path = os.path.join("figs/ICE_PIC/data/case_CH4_counterflow_premixed", target_case, "predicted_Y.csv")
ice_data_path = os.path.join("figs/GFRI/data/PaSR", target_case, "predicted_Y.csv")

# Read data
df_original = pd.read_csv(original_data_path)
df_rcce = pd.read_csv(rcce_data_path)
df_gfri = pd.read_csv(gfri_data_path)
df_ice = pd.read_csv(ice_data_path)

MF = df_original['MF']  # Adjust 'x' column name if necessary

################### Set figure and axis layout
fig, axs = plt.subplots(3, 3, figsize=(15/inch, 12/inch), sharex=True)

# Define columns for OH, CH2O, HCO, CEM, MF, HRR (Qdot)
data_cols = ['OH', 'O','H','CH2O', 'HCO','CH','CEM', 'MF', 'Qdot']
data_list = [(df_original[col], df_rcce[col], df_gfri[col], df_ice[col], col) for col in data_cols]

# Plot data
for i, (original, rcce, gfri, ice, label) in enumerate(data_list):
    row, col = divmod(i, 3)
    ax = axs[row, col]

    # Scatter plot for original data
    ax.scatter(MF, original, label='Original', color='black', s=3, marker='o', alpha=0.5)
    
    # Scatter plot for RCCE, GFRI, and ICE
    # ax.scatter(MF, rcce, label='RCCE', color='blue', s=3, marker='x', alpha=0.3)
    ax.scatter(MF, gfri, label='GFRI', color='red', s=3, marker='*', alpha=0.3)
    # ax.scatter(MF, ice, label='ICE', color='green', s=10, marker='*', alpha=0.3)
    
    # Add labels and limits
    ax.set_ylabel(label, fontsize=10)
    # ax.set_xscale('log')  # Set x-axis scale to logarithmic
    
    if row == 2:
        ax.set_xlabel('x', fontsize=10)

# Hide any unused subplots
if len(data_list) < 9:
    for j in range(len(data_list), 9):
        fig.delaxes(axs[j // 3, j % 3])

# Apply formatting and axis adjustments
for ax_row in axs:
    for ax in ax_row:
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        # ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))  
        ax.xaxis.get_major_formatter().set_powerlimits((-2, 2))
        ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
        ax.set_xlim(x_range)
        
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.92, 0.92), fontsize=8, borderaxespad=0)
plt.tight_layout()
plt.subplots_adjust(wspace=0.5, hspace=0.2)
fig_dir = "figs/combined_plot"
os.makedirs(fig_dir, exist_ok=True)
plt.savefig(os.path.join(fig_dir, "PaSR_CEQ_combine.png"), dpi=300)

################### figure2 comparison gfri #####################
############# figure2 comparison gfri #####################
# Define the data columns to compare
data_cols = ['OH', 'O', 'H', 'CH2O', 'HCO', 'CH', 'CEM', 'MF', 'Qdot']

# Create a new figure with 9 subplots (3x3 grid)
fig, axs = plt.subplots(3, 3, figsize=(8, 8))

# Flatten the axs array for easy iteration
axs = axs.flatten()

# Loop through each data column and create a subplot
for i, col in enumerate(data_cols):
    ax = axs[i]
    
    # Scatter plot of original[col] vs gfri[col], colored by df_original['T']
    scatter = ax.scatter(df_original[col], df_gfri[col], c=df_original['T'], cmap='jet', s=2, alpha=0.8)
    
    # Plot y=x line
    min_val = min(df_original[col].min(), df_gfri[col].min())
    max_val = max(df_original[col].max(), df_gfri[col].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')
    
    # Add labels and title
    ax.set_title(f'{col}')
    
    # Set axis limits to be equal
    ax.set_xlim([min_val, max_val])
    ax.set_ylim([min_val, max_val])
    
    # Set major formatter for x and y axes
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
    ax.xaxis.get_major_formatter().set_powerlimits((-2, 2))
    ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))       

# cbar = fig.colorbar(scatter, ax=axs, orientation='vertical', fraction=0.02, pad=0.04)
# cbar.set_label('Temperature (T)')

plt.tight_layout()             

# Save the figure
fig_dir = "figs/combined_plot"
os.makedirs(fig_dir, exist_ok=True)
plt.savefig(os.path.join(fig_dir, "Original_vs_GFRI_comparison.png"), dpi=300)

# ################### figure3

# # Create a new figure with 9 subplots (3x3 grid)
# fig, axs = plt.subplots(3, 3, figsize=(8, 8))

# # Flatten the axs array for easy iteration
# axs = axs.flatten()

# # Loop through each data column and create a subplot
# for i, col in enumerate(data_cols):
#     ax = axs[i]
    
#     # Scatter plot of original[col] vs rcce[col]
#     ax.scatter(df_original[col], df_rcce[col], label=f'Original vs RCCE {col}', c=df_original['T'], cmap='jet', s=2, alpha=0.8)
    
#     # Plot y=x line
#     min_val = min(df_original[col].min(), df_rcce[col].min())
#     max_val = max(df_original[col].max(), df_rcce[col].max())
#     ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')
    
#     # Add labels and title
#     # ax.set_xlabel(f'Original {col}')
#     # ax.set_ylabel(f'RCCE {col}')
#     ax.set_title(f'{col}')
#     # ax.legend()
    
#     # Set axis limits to be equal
#     ax.set_xlim([min_val, max_val])
#     ax.set_ylim([min_val, max_val])
    
#     # Set major formatter for x and y axes
#     ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#     ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#     ax.xaxis.get_major_formatter().set_powerlimits((-2, 2))
#     ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))

# # Adjust layout to prevent overlap
# plt.tight_layout()

# # Save the figure
# fig_dir = "figs/combined_plot"
# os.makedirs(fig_dir, exist_ok=True)
# plt.savefig(os.path.join(fig_dir, "Original_vs_RCCE_comparison.png"), dpi=300)

