import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

inch = 2.54
target_case = "PSR_res"
original_data_path = os.path.join("data/PSR", target_case) + ".csv"
rcce_data_path = os.path.join("figs/RCCE/data/PSR", target_case, "predicted_Y.csv")
gfri_data_path = os.path.join("figs/GFRI/data/PSR", target_case, "predicted_Y.csv")
# ice_data_path = os.path.join("figs/ICE_PIC/data/case_CH4_counterflow_premixed", target_case, "predicted_Y.csv")
ice_data_path = os.path.join("figs/GFRI/data/PSR", target_case, "predicted_Y.csv")



# Read data
df_original = pd.read_csv(original_data_path)
df_rcce = pd.read_csv(rcce_data_path)
df_gfri = pd.read_csv(gfri_data_path)
df_ice = pd.read_csv(ice_data_path)

t_res = df_original['t_res']  # Adjust 'x' column name if necessary

# Set figure and axis layout
fig, axs = plt.subplots(3, 3, figsize=(15/inch, 12/inch), sharex=True)

# Define columns for OH, CH2O, HCO, CEM, MF, HRR (Qdot)
data_cols = ['OH', 'O','H','CH2O', 'HCO','CH','CEM', 'MF', 'Qdot']
data_list = [(df_original[col], df_rcce[col], df_gfri[col], df_ice[col], col) for col in data_cols]

# Plot data
for i, (original, rcce, gfri, ice, label) in enumerate(data_list):
    row, col = divmod(i, 3)
    ax = axs[row, col]

    # Scatter plot for original data
    ax.scatter(t_res, original, label='Original', color='blue', s=10, marker='x', alpha=0.8)
    
    # Line plot for RCCE, GFRI, and ICE
    ax.plot(t_res, rcce, label='RCCE', linestyle='--', color='black', linewidth=1)
    ax.plot(t_res, gfri, label='GFRI', linestyle=':', color='red', linewidth=1)
    ax.plot(t_res, ice, label='ICE', linestyle='-', color='green', linewidth=1)
    
    # Add labels and limits
    ax.set_ylabel(label, fontsize=10)
    ax.set_xlim([t_res.min(), t_res.max()])
    ax.set_xscale('log')  # Set x-axis scale to logarithmic
    
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
handles, labels = axs[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.92, 0.92), fontsize=8, borderaxespad=0)
plt.tight_layout()
plt.subplots_adjust(wspace=0.5, hspace=0.2)
fig_dir = "figs/combined_plot"
os.makedirs(fig_dir, exist_ok=True)
plt.savefig(os.path.join(fig_dir, "PSR_CEQ_combine.png"), dpi=300)

