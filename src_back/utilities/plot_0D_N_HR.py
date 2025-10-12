import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D

QSSA_path = "figs/QSSA/data/"
RCCE_path = "figs/RCCE/data/"

case_name = "case_NH3_0D"

line_arr = ('-', '--', '-.', ':')
color_arr = ('k', 'r', 'b', 'y', 'g', 'c', 'm')
symbol_arr = ('s', 'o', 'v', '^', '*')

# T_range = [1500,3600]

########## FIGURE 1: PLOT against Time ###########
def plot_over_time_compact_ori():
    # Define target cases
    QSSA_path = "figs/QSSA/data/"
    RCCE_path = "figs/RCCE/data/"
    case_name = "case_NH3_0D"

    line_arr = ('-', '--', '-.', ':')
    color_arr = ('k', 'r', 'b', 'y', 'g', 'c', 'm')
    symbol_arr = ('s', 'o', 'v', '^', '*')

    target_case_assemble = [
        'N_HR_1',
        'N_HR_2',
        'N_HR_3',
    ]

    # Corresponding phi labels for legend
    name_case_assemble = [
        r'$\phi = 0.8$',
        r'$\phi = 1.0$',
        r'$\phi = 1.2$',
    ]

    # Define data columns for the top row
    top_row_left_cols = ['NH3', 'O2', 'H2', 'H2O']
    top_row_right_col = 'T'

    # Define data columns for the lower rows
    data_cols = ['OH', 'NO', 'CEM', 'Qdot']

    x_range = [1e-5, 1e-2]
    range_cols = {'Qdot': (1e6, 1e12)}
    log_cols = {'Qdot': True}

    plt.rcParams.update({'font.size': 8})

    # Create subplots: len(data_cols) + 1 rows and len(target_case_assemble) columns
    fig, axs = plt.subplots(len(data_cols) + 1, len(target_case_assemble),
                            figsize=(4.73, 1.3 * (len(data_cols) + 1))
                            )

    # Create lists to store legend handles and labels for the first row
    first_row_handles = []
    first_row_labels = []

    # Loop through cases for the top row
    for col_idx, target_case in enumerate(target_case_assemble):
        # File paths for the current case
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"

        # Read data
        df_original = pd.read_csv(original_data_path)
        t = df_original['t'].to_numpy()

        # Access the top-row subplot
        ax = axs[0, col_idx]
        ax_right = ax.twinx()

        # Plot species on the left y-axis
        for col in top_row_left_cols:
            line, = ax.plot(
                t[::1],
                df_original[col].to_numpy()[::1],
                label=col,
                alpha=1,
                linestyle='-',
                linewidth=1.0
            )
            if col_idx == 0:  # Add legend for only the first column
                first_row_handles.append(line)
                first_row_labels.append(col)

        # Plot temperature on the right y-axis
        line_temp, = ax_right.plot(
            t[::1],
            df_original[top_row_right_col].to_numpy()[::1],
            label=top_row_right_col,
            alpha=1,
            linestyle='--',
            linewidth=1.0,
            color='black'
        )
        if col_idx == 0:  # Add legend for temperature for only the first column
            first_row_handles.append(line_temp)
            first_row_labels.append(top_row_right_col)

        # Format axes
        ax.set_xscale('log')
        ax.set_xlim(x_range)
        ax.set_xticklabels([])  # To remove y-axis labels for non-bottom rows
        ax.set_title(name_case_assemble[col_idx], fontsize=8)
        
        # Add "T" as the title for the rightmost subplot
        if col_idx == 0:
            ax.set_ylabel("X", fontsize=8)


        if col_idx == len(target_case_assemble) - 1:
            ax_right.set_ylabel("T (K)", fontsize=8)
            # ax_right.set_title("T", fontsize=8, loc='right')

        # Set xticks and yticks only for the leftmost and rightmost subplots
        if col_idx != 0:
            ax.set_yticks([])
        if col_idx != len(target_case_assemble) - 1:
            ax_right.set_yticks([])

        # Set x-axis ticks for only the bottom row subplots
        if col_idx == len(target_case_assemble) - 1:
            ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2])

    # Add legend for all subplots in the first row
    for col_idx in range(len(target_case_assemble)):
        ax = axs[0, col_idx]
        legend_1 = ax.legend(
            first_row_handles,
            first_row_labels,
            loc='upper center',
            ncol=1,
            fontsize=6,
            frameon=False,
            bbox_to_anchor=(0.8, 0.95)  # Adjust the positioning to avoid overlap
        )
        # Synchronize color of legend texts with line colors
        for text, line in zip(legend_1.get_texts(), legend_1.get_lines()):
            text.set_color(line.get_color())

    # Loop through cases for the lower rows
    for col_idx, target_case in enumerate(target_case_assemble):
        # File paths for the current case
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"
        rcce_data_path = os.path.join(RCCE_path, case_name, target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(QSSA_path, case_name, target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        t = df_original['t'].to_numpy()

        for row_idx, col in enumerate(data_cols):
            ax = axs[row_idx + 1, col_idx]  # Offset by 1 for the top row

            # Plot CVODE (line plot)
            # --- CVODE (the "true" data) ---
            cvode_values = df_original[col].to_numpy()
            ax.plot(
                t,
                cvode_values,
                linestyle='-',       # solid line
                linewidth=1,       # slightly thicker
                color='black',       # black line
                label='CVODE' if col_idx == 0 and row_idx == 0 else ""
            )

            # --- QSSA (predicted) as red, open circles ---
            qssa_values = df_qssa[col].to_numpy()
            ax.scatter(
                t[::50],
                qssa_values[::50],
                s=5,                # larger markers than default
                alpha=0.8,
                marker='o',          # circle
                edgecolors='#D55E00',   # orange outline
                facecolors='none',      # open circle
                # edgecolors='purple',    # red edges
                label='QSSA' if col_idx == 0 and row_idx == 0 else ""
            )

            # --- RCCE (predicted) as blue squares ---
            rcce_values = df_rcce[col].to_numpy()
            ax.scatter(
                t[::50],
                rcce_values[::50],
                s=5,
                alpha=0.8,
                marker='s',          # square
                edgecolors='#0072B2',   # sky blue outline
                # edgecolors='#56B4E9',   # sky blue outline
                facecolors='none',      # open circle
                # color='green',        # solid blue squares
                label='RCCE' if col_idx == 0 and row_idx == 0 else ""
            )

            # Format axes
            ax.set_xscale('log')
            ax.set_xlim(x_range)
            
            if col == 'Qdot':
                ax.set_yscale('log')
            
            if row_idx + 1 == len(data_cols):  # Only bottom row gets x-axis label
                ax.set_xlabel("t (s)", fontsize=8)
            if col_idx == 0:
                if col == 'Qdot':
                    ax.set_ylabel('HRR (J/kg/s)', fontsize=8)
                else:
                    ax.set_ylabel(col, fontsize=8)

            # Set y-axis ticks for non-first columns
            if col_idx != 0:
                # ax.set_yticks([])
                ax.set_yticklabels([])  # To remove y-axis labels for non-bottom rows
            
            # Remove x-axis labels for non-bottom rows
            if row_idx != len(data_cols) - 1:
                ax.set_xticklabels([])

    legend_handles = [
        Line2D([0], [0],
            linestyle=line_arr[0],
            linewidth=1.0,
            color='black',
            label='CVODE'
        ),
        Line2D([0], [0],
            marker='o',
            linestyle='None',
            markerfacecolor='none',   # use markerfacecolor, not facecolors
            markeredgecolor='#D55E00',
            alpha=1,
            label='QSSA'
        ),
        Line2D([0], [0],
            marker='s',
            linestyle='None',
            markerfacecolor='none',
            markeredgecolor='#0072B2',
            alpha=1,
            label='CEQ'
        )
    ]


    plt.subplots_adjust(top=0.90, right=0.85, wspace=0.1, hspace=0.1)

    # Add the legend on top, centered horizontally
    fig.legend(
        handles=legend_handles,
        loc='upper center',
        ncol=3,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, 0.98)
    )

    # Save the figure
    fig.savefig("figs/0D_compact.png", dpi=300, bbox_inches='tight')
    fig.savefig("figs/0D_compact.pdf", dpi=300, bbox_inches='tight')
    plt.show()



########## FIGURE 2: PLOT against Temperature ###########
def plot_over_temperature():
        
    # Define target cases
    target_case_assemble = [
        # 'N_HR_1',
        'N_HR_2',
        # 'N_HR_3',
        # 'N_HR_4',
        # 'N_HR_5',
        # 'N_HR_6',
        # 'N_HR_7',
        # 'N_HR_8',
        # 'N_HR_9'
    ]

    name_case_assemble = [
        # r'$\phi = 0.8$',
        r'$\phi = 1.0$',
        # r'$\phi = 1.2$',
    ]

    # Define data columns to plot
    data_cols = ['OH', 'NO', 'CEM', 'Qdot']

    # Define y-axis ranges for each data column (optional)
    range_cols = {
        'Qdot': (1e6, 1e12)  # Example range for Qdot
    }

    # Define which columns should use a logarithmic y-axis
    log_cols = {
        'Qdot': True  # Use log scale for Qdot
    }

    # Create subplots: rows = data_cols, columns = target_case_assemble
    fig, axs = plt.subplots(len(data_cols), len(target_case_assemble), figsize=(5.67, 7), sharex=False)

    # Ensure axs is a 2D array even if there's only one row or column
    if len(data_cols) == 1 and len(target_case_assemble) == 1:
        axs = np.array([[axs]])
    elif len(data_cols) == 1:
        axs = axs[np.newaxis, :]
    elif len(target_case_assemble) == 1:
        axs = axs[:, np.newaxis]

    # Store x-range for each target_case to ensure consistency
    x_ranges = {}

    # Iterate over each target case
    for j, target_case in enumerate(target_case_assemble):
        # Define file paths
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"
        rcce_data_path = os.path.join(RCCE_path, case_name, target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(QSSA_path, case_name, target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        # Extract temperature data
        T = df_original['T'].to_numpy()

        # Store the x-range (temperature range) for the current target_case
        x_range = (min(T), max(T))
        x_ranges[target_case] = x_range

        # Iterate over each data column
        for i, col in enumerate(data_cols):
            ax = axs[i, j]  # Access the correct subplot

            # Scatter plot for original data ("VODE")
            scatter_plot = ax.scatter(
                T[::60],
                df_original[col].to_numpy()[::60],
                s=4,
                alpha=0.8,
                marker=symbol_arr[0],
                color=color_arr[0],
                label='VODE' if (i == 0 and j == 0) else ""
            )

            # Line plot for QSSA data
            qssa_plot = ax.plot(
                T,
                df_qssa[col].to_numpy(),
                linestyle=line_arr[1],
                linewidth=1.0,
                color=color_arr[1],
                label='QSSA' if (i == 0 and j == 0) else ""
            )

            # Line plot for RCCE data
            rcce_plot = ax.plot(
                T,
                df_rcce[col].to_numpy(),
                linestyle=line_arr[2],
                linewidth=1.0,
                color=color_arr[2],
                label='RCCE' if (i == 0 and j == 0) else ""
            )

            # Set logarithmic y-axis if specified
            if log_cols.get(col, False):
                ax.set_yscale('log')
                # Set fixed tick positions at the desired powers of 10
                # ax.yaxis.set_major_locator(mticker.FixedLocator([1e6, 1e8, 1e10, 1e12]))
                # Use a log formatter to show them as 10^6, 10^8, etc.
                ax.yaxis.set_major_formatter(mticker.LogFormatterSciNotation(base=10))
            else: 
                ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                # ax.yaxis.get_major_formatter().set_powerlimits((-3, 3))
                
            if col in range_cols:
                ax.set_ylim(range_cols[col])

            # Set x-axis limits consistent for the same target_case
            ax.set_xlim(x_ranges[target_case])

            # Set x-axis ticks to appear every 500 units
            ax.xaxis.set_major_locator(mticker.MultipleLocator(500))

            # Set y-label only on the first column
            if j == 0:
                if col == 'Qdot':
                    ax.set_ylabel('HRR (J/kg/s)', fontsize=10)
                else:
                    ax.set_ylabel(col, fontsize=10)
            else:
                ax.set_ylabel('')  # Remove y-label for non-first columns



            # Set x-label only on the bottom row
            if i == len(data_cols) - 1:
                ax.set_xlabel('T(K)', fontsize=10)

        # Add title for each target_case at the top of each column
        axs[0, j].set_title(name_case_assemble[j], fontsize=10)

    # Create custom legend handles
    legend_handles = [
        Line2D([0], [0], marker=symbol_arr[0], color='w', label='VODE',
            markerfacecolor=color_arr[0], markersize=6, linestyle='None'),
        Line2D([0], [0], color=color_arr[1], linestyle=line_arr[1], linewidth=1.0, label='QSSA'),
        Line2D([0], [0], color=color_arr[2], linestyle=line_arr[2], linewidth=1.0, label='RCCE')
    ]

    plt.subplots_adjust(top=0.90, right=0.85, wspace=0.3, hspace=0.4)

    # Add the legend on top, centered horizontally
    fig.legend(
        handles=legend_handles,
        labels=['VODE', 'QSSA', 'RCCE'],
        loc='upper center',
        ncol=3,
        fontsize=10,
        title_fontsize='11',
        frameon=False
    )

    # Maintain tight layout while accommodating the legend
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # Save the figure
    fig_dir = "figs/NH3_0Dignition"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "QSSA_RCCE_VODE_vs_Original_Temperature_PerCase.png"), dpi=300)
    plt.show()



########## FIGURE 3: PLOT against Time (Compact) ###########

def plot_over_time_compact():
    # Define target cases
    target_case_assemble = [
        'N_HR_1',
        'N_HR_2',
        'N_HR_3',
    ]

    # Define names for cases (used for column titles)
    name_case_assemble = [
        r'$\phi = 0.8$',
        r'$\phi = 1.0$',
        r'$\phi = 1.2$'
    ]

    # Define data columns to plot
    data_cols = ['OH', 'NO', 'CEM', 'Qdot']

    x_range = [1e-5, 1e-2]
    range_cols = {'Qdot': (1e6, 1e12)}
    log_cols = {'Qdot': True}

    plt.rcParams.update({'font.size': 8})

    # Create a grid of subplots with rows=len(data_cols) and cols=len(target_case_assemble)
    fig, axs = plt.subplots(len(data_cols), len(target_case_assemble),
                            figsize=(1.5 * len(target_case_assemble), 1.5 * len(data_cols)),
                            sharex=True, sharey='row')

    # Loop through each case and plot
    for col_idx, target_case in enumerate(target_case_assemble):
        # File paths for the current case
        original_data_path = os.path.join("data", case_name,target_case) + ".csv"
        rcce_data_path = os.path.join(RCCE_path, case_name,target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(QSSA_path, case_name, target_case,"predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        t = df_original['t'].to_numpy()

        for row_idx, col in enumerate(data_cols):
            ax = axs[row_idx, col_idx]  # Access subplot at (row_idx, col_idx)

            # Plot CVODE (line plot)
            cvode_values = df_original[col].to_numpy()
            ax.plot(
                t,
                cvode_values,
                linestyle=line_arr[0],
                linewidth=1.0,
                color=color_arr[0],
                label='CVODE' if col_idx == 0 and row_idx == 0 else ""
            )

            # Find max value and position for CVODE
            cvode_max_idx = np.argmax(cvode_values)
            cvode_max_value = cvode_values[cvode_max_idx]
            cvode_max_time = t[cvode_max_idx]

            # Plot QSSA (scatter)
            qssa_values = df_qssa[col].to_numpy()
            ax.scatter(
                t[::60],
                qssa_values[::60],
                s=1,
                alpha=0.8,
                marker=symbol_arr[1],
                color=color_arr[1],
                label='QSSA' if col_idx == 0 and row_idx == 0 else ""
            )

            # Plot RCCE (scatter)
            rcce_values = df_rcce[col].to_numpy()
            ax.scatter(
                t[::60],
                rcce_values[::60],
                s=1,
                alpha=0.8,
                marker=symbol_arr[2],
                color=color_arr[2],
                label='RCCE' if col_idx == 0 and row_idx == 0 else ""
            )

            # Set logarithmic y-axis if specified
            if log_cols.get(col, False):
                ax.set_yscale('log')
                ax.yaxis.set_major_formatter(mticker.LogFormatterSciNotation(base=10))
            else:
                ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))

            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(mticker.LogFormatterSciNotation(base=10))
            ax.set_xlim(x_range)

            if col in range_cols:
                ax.set_ylim(range_cols[col])

            # Set y-axis labels for the first column only
            if col_idx == 0:
                ax.set_ylabel('HRR (J/kg/s)' if col == 'Qdot' else col, fontsize=6)

            # Set x-axis labels for the last row only
            if row_idx == len(data_cols) - 1:
                ax.set_xlabel('t (s)', fontsize=6)

            # Add column titles
            if row_idx == 0:
                ax.set_title(name_case_assemble[col_idx], fontsize=7, pad=10)

    # Create custom legend
    legend_handles = [
        Line2D([0], [0], color=color_arr[0], linestyle=line_arr[0], linewidth=1.0, label='CVODE'),
        Line2D([0], [0], marker=symbol_arr[1], color='w', label='QSSA',
               markerfacecolor=color_arr[1], markersize=4, linestyle='None'),
        Line2D([0], [0], marker=symbol_arr[2], color='w', label='RCCE',
               markerfacecolor=color_arr[2], markersize=4, linestyle='None')
    ]

    fig.legend(
        handles=legend_handles,
        labels=['VODE', 'QSSA', 'RCCE'],
        loc='upper center',
        ncol=3,
        fontsize=8,
        title_fontsize=8,
        frameon=False
    )

    plt.tight_layout(rect=[0, 0, 1, 0.9])

    # Save the figure
    fig_dir = "figs/NH3_0Dignition"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "QSSA_RCCE_VODE_vs_Original_Time_MultiCase.png"), dpi=300)
    plt.show()




###############MAIN##################
def main():
    data_cols = ['OH', 'O', 'H', 'N', 'NO', 'NO2', 'CEM', 'MF', 'Qdot']
    plot_over_time_compact_ori()
    plot_over_time_compact()
    


if __name__ == "__main__":
    main()
    