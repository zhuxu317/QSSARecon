import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D

qssa_method_path = "figs/QSSA/data/"
# rcce_method_path = "figs/RCCE/data/"
rcce_method_path = "figs/ILDM/data/"

case_name = "case_NH3_counterflow"
line_arr = ('-','--','-.',':')
color_arr = ('k','r','b','y','g', 'c','m')
symbol_arr = ('s','o','v','^','*')



###########Figure 1: compare the QSSA and RCCE over T overall##############
# Define target cases
def plot_results_over_T():
    target_case_assemble = [
        'N_CF_1',
        # 'N_CF_2',
        # 'N_CF_3',
        # 'N_CF_4',
        # 'N_CF_5',
        # 'N_CF_6',
        # 'N_CF_7',
        # 'N_CF_8',
    ]

    x_range_assemble = [
        [18, 27],
        # [18, 35],
        # [17, 26],
        # [17, 26],
        # [18, 35],
        # [22, 31],
        # [22, 31],
        # [22, 31],
    ]

    # Initialize subplots
    fig, axs = plt.subplots(3, 3, figsize=(8, 7), sharex=True)

    data_cols = ['OH', 'O', 'H', 'N', 'NO', 'NO2', 'CEM', 'MF', 'Qdot']

    # Define color and line styles for QSSA and RCCE
    qssa_color = color_arr[1] #'blue'
    rcce_color = color_arr[2]# 'red'
    qssa_linestyle = '-'
    rcce_linestyle = '--'

    for j, target_case in enumerate(target_case_assemble):
        # Define paths
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"
        rcce_data_path = os.path.join(rcce_method_path, case_name, target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(qssa_method_path, case_name, target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        # Ensure T is a NumPy array
        T = df_original['T'].to_numpy()

        # Plot each variable in its corresponding subplot
        for i, col in enumerate(data_cols):
            row, col_num = divmod(i, 3)
            ax = axs[row, col_num]
            
            # Scatter plot for original data
            ax.scatter(T[::3], 
                    df_original[col].to_numpy()[::3],  # Convert to NumPy array
                    label=f'{target_case}', 
                    s=10, 
                    alpha=0.8, 
                    marker='o', 
                    color='black')  # Original data in black
            
            # Line plot for QSSA
            ax.plot(T, 
                    df_qssa[col].to_numpy(),  # Convert to NumPy array
                    linestyle=qssa_linestyle, 
                    linewidth=1.5, 
                    color=qssa_color, 
                    label=f'QSSA - {target_case}')
            
            # Line plot for RCCE
            ax.plot(T, 
                    df_rcce[col].to_numpy(),  # Convert to NumPy array
                    linestyle=rcce_linestyle, 
                    linewidth=1.5, 
                    color=rcce_color, 
                    label=f'RCCE - {target_case}')

            # Formatting
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
            ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
            ax.set_ylabel(col, fontsize=10)
            ax.set_xlim([1000, 2500])
            if row == 2:
                ax.set_xlabel('T(K)', fontsize=10)

    # Remove extra subplots
    for j in range(len(data_cols), 9):
        fig.delaxes(axs[j // 3, j % 3])

    # # Add a legend
    # fig.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fontsize=10)

    plt.tight_layout()

    # Save the figure
    fig_dir = "figs/NH3"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_Temperature.png"), dpi=300)
    plt.show()


###########Figure 2: compare the QSSA and RCCE over T overall##############
# Define target cases
def plot_results_over_T_overall():
    target_case_assemble = [
        'N_CF_1',
        # 'N_CF_2',
        # 'N_CF_3',
        # 'N_CF_4',
        # 'N_CF_5',
        # 'N_CF_6',
        # 'N_CF_7',
        # 'N_CF_8',
    ]

    x_range_assemble = [
        [18, 27],
        # [18, 35],
        # [17, 26],
        # [17, 26],
        # [18, 35],
        # [22, 31],
        # [22, 31],
        # [22, 31],
    ]

    # Initialize subplots (2 rows and 3 columns)
    fig, axs = plt.subplots(2, 3, figsize=(5.73, 3), sharex=True)

    data_cols = ['MF', 'Qdot', 'CEM', 'OH', 'NO2', 'NO']

    # Define color and line styles for QSSA and RCCE
    qssa_color = color_arr[1]  # 'blue'
    rcce_color = color_arr[2]  # 'red'
    qssa_linestyle = '-'
    rcce_linestyle = '--'

    # Define y-axis ranges for each data column (optional)
    range_cols = {
        'Qdot': (1e6, 1e12),  # Example range for Qdot (HRR)
    }

    for j, target_case in enumerate(target_case_assemble):
        # Define paths
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"
        rcce_data_path = os.path.join(rcce_method_path, case_name, target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(qssa_method_path, case_name, target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        # Ensure T is a NumPy array
        T = df_original['T'].to_numpy()

        # Plot each variable in its corresponding subplot
        for i, col in enumerate(data_cols):
            row, col_num = divmod(i, 3)
            ax = axs[row, col_num]
            
            # Scatter plot for original data
            ax.scatter(T[::3], 
                    df_original[col].to_numpy()[::3],  # Convert to NumPy array
                    label=f'{target_case}', 
                    s=10, 
                    alpha=0.8, 
                    marker='o', 
                    color='black')  # Original data in black
            
            # Line plot for QSSA
            ax.plot(T, 
                    df_qssa[col].to_numpy(),  # Convert to NumPy array
                    linestyle=qssa_linestyle, 
                    linewidth=1.5, 
                    color=qssa_color, 
                    label=f'QSSA - {target_case}')
            
            # Line plot for RCCE
            ax.plot(T, 
                    df_rcce[col].to_numpy(),  # Convert to NumPy array
                    linestyle=rcce_linestyle, 
                    linewidth=1.5, 
                    color=rcce_color, 
                    label=f'RCCE - {target_case}')

            # Log scale for Qdot (HRR) only
            if col == 'Qdot':
                ax.set_yscale('log')
                ax.set_ylabel('HRR (J/kg/s)', fontsize=10)  # Update label for Qdot
                if col in range_cols:  # Apply the custom range if specified
                    ax.set_ylim(range_cols[col])

            # Formatting for other axes
            else:
                ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
                ax.set_ylabel(col, fontsize=10)

            ax.set_xlim([1000, 2300])
            if row == 1:  # Add x-axis label only on the bottom row
                ax.set_xlabel('T(K)', fontsize=10)

    # Custom legend
    legend_handles = [
        Line2D([0], [0], marker='s', color='w', label='Original Data',
            markerfacecolor='black', markersize=4, linestyle='None'),
        Line2D([0], [0], color=qssa_color, linestyle=qssa_linestyle, linewidth=1.5, label='QSSA'),
        Line2D([0], [0], color=rcce_color, linestyle=rcce_linestyle, linewidth=1.5, label='RCCE'),
    ]

    fig.legend(
        handles=legend_handles, 
        labels=['VODE','QSSA', 'RCCE'],
        loc='upper center', 
        bbox_to_anchor=(0.5, 0.99), 
        ncol=3, 
        frameon=False)

    # fig.tight_layout()
    plt.subplots_adjust(wspace=0.65, hspace=0.2, top=0.89, bottom=0.15)

    # Save the figure
    fig_dir = "figs/NH3"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_Temperature_overall.png"), dpi=300)
    plt.show()


# ###################### Plot FIGURE RECONSTRUCT SPACE ######################
# Define target cases
def plot_results_over_space():

    target_case_assemble = [
        'N_CF_1',
        # 'N_CF_2',
        # 'N_CF_3',
        # 'N_CF_4',
        # 'N_CF_5',
        # 'N_CF_6',
        # 'N_CF_7',
        # 'N_CF_8',
    ]

    x_range_assemble = [
        [17, 28],
        # [18, 35],
        # [17, 26],
        # [17, 26],
        # [18, 35],
        # [22, 31],
        # [22, 31],
        # [22, 31],
    ]

    # Initialize subplots (2 rows and 3 columns)
    fig, axs = plt.subplots(2, 3, figsize=(4.73, 3), sharex=True)

    data_cols = ['MF', 'Qdot', 'CEM', 'OH', 'NO2', 'NO']

    # Define color and line styles for QSSA and RCCE
    qssa_color = color_arr[1]  # 'blue'
    rcce_color = color_arr[2]  # 'red'
    qssa_linestyle = '-'
    rcce_linestyle = '--'

    # Define y-axis ranges for each data column (optional)
    range_cols = {
        'Qdot': (1e6, 1e12),  # Example range for Qdot (HRR)
    }

    for j, target_case in enumerate(target_case_assemble):
        # Define paths
        original_data_path = os.path.join("data", case_name, target_case) + ".csv"
        rcce_data_path = os.path.join(rcce_method_path, case_name, target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(qssa_method_path, case_name, target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        # Ensure x is a NumPy array (Changed from 'T' to 'x')
        x = df_original['grid'].to_numpy() * 1000  # Update column name if different

        # Plot each variable in its corresponding subplot
        for i, col in enumerate(data_cols):
            row, col_num = divmod(i, 3)
            ax = axs[row, col_num]

            # Scatter plot for original data
            ax.scatter(
                x[::3],  # Changed from T to x
                df_original[col].to_numpy()[::3],  # Convert to NumPy array
                label=f'{target_case}',
                s=10,
                alpha=0.8,
                marker='o',
                color='black'  # Original data in black
            )

            # Line plot for QSSA
            ax.plot(
                x,  # Changed from T to x
                df_qssa[col].to_numpy(),  # Convert to NumPy array
                linestyle=qssa_linestyle,
                linewidth=1.5,
                color=qssa_color,
                label=f'QSSA - {target_case}'
            )

            # Line plot for RCCE
            ax.plot(
                x,  # Changed from T to x
                df_rcce[col].to_numpy(),  # Convert to NumPy array
                linestyle=rcce_linestyle,
                linewidth=1.5,
                color=rcce_color,
                label=f'RCCE - {target_case}'
            )

            # Log scale for Qdot (HRR) only
            if col == 'Qdot':
                ax.set_yscale('log')
                ax.set_ylabel('HRR (J/kg/s)', fontsize=10)  # Update label for Qdot
                if col in range_cols:  # Apply the custom range if specified
                    ax.set_ylim(range_cols[col])

            # Formatting for other axes
            else:
                ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
                ax.set_ylabel(col, fontsize=10)

            # Set x-axis limits using x_range_assemble
            ax.set_xlim(x_range_assemble[j])  # Changed from [1000, 2300] to x range

            if row == 1:  # Add x-axis label only on the bottom row
                ax.set_xlabel('x (mm)', fontsize=10)  # Updated label from 'T(K)' to 'x (units)'

    # Custom legend
    legend_handles = [
        Line2D([0], [0], marker='s', color='w', label='VODE',
            markerfacecolor='black', markersize=4, linestyle='None'),
        Line2D([0], [0], color=qssa_color, linestyle=qssa_linestyle, linewidth=1.5, label='QSSA'),
        Line2D([0], [0], color=rcce_color, linestyle=rcce_linestyle, linewidth=1.5, label='RCCE'),
    ]

    fig.legend(
        handles=legend_handles,
        labels=['VODE', 'QSSA', 'RCCE'],  # Corrected label from 'VODE' to 'Original Data'
        loc='upper center',
        bbox_to_anchor=(0.5, 0.99),
        ncol=3,
        frameon=False
    )

    # Adjust layout
    plt.subplots_adjust(wspace=0, hspace=0.2, top=0.89, bottom=0.15)
    # Save the figure
    fig_dir = "figs/NH3"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_x_overall.png"), dpi=300)  # Updated filename
    plt.show()
    
def plot_results_over_space_noise():
    target_case = 'N_CF_1'
    target_noise_assemble = [0.0, 0.15, 0.3, 0.45]
    x_range = [20.5, 25]

    # Initialize subplots (len(data_cols) rows, len(target_noise_assemble) columns)
    data_cols = ['MF', 'CEM', 'OH', 'NO']

    fig, axs = plt.subplots(
        len(data_cols), 
        len(target_noise_assemble), 
        figsize=(5.7, 5), 
        sharex=True
    )

    # Define colors for QSSA and RCCE
    qssa_color = 'red'
    rcce_color = 'steelblue'
    # rcce_color = 'blue'
    

    # Define y-axis ranges for each data column (optional)
    range_cols = {}


    # Loop through target_noise_assemble values
    for j, target_noise in enumerate(target_noise_assemble):
        # Define paths
        base_case_name = 'case_NH3_counterflow'
        original_data_path = os.path.join("data", base_case_name, target_case) + ".csv"
        noise_target_case = f'{target_case}_{target_noise}'

        noise_case_name = 'case_NH3_counterflow_noise'
        rcce_data_path = os.path.join(rcce_method_path, noise_case_name, noise_target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(qssa_method_path, noise_case_name, noise_target_case, "predicted_X.csv")

        # Read data
        df_original = pd.read_csv(original_data_path)
        df_rcce = pd.read_csv(rcce_data_path)
        df_qssa = pd.read_csv(qssa_data_path)

        # Ensure x is a NumPy array
        x = df_original['grid'].to_numpy() * 1000  # Convert to appropriate units if necessary

        # Plot each variable in its corresponding subplot
        for i, col in enumerate(data_cols):
            ax = axs[i, j]  # rows = data_cols, columns = target_noise

            # Scatter plot for Original Data (VODE)
            ax.plot(
                x[::2],  # Sample every 3rd point for clarity
                df_original[col].to_numpy()[::2],
                label='VODE',
                alpha=1,
                color='black'
            )

            # Scatter plot for QSSA with hollow markers
            ax.scatter(
                x[::1],
                df_qssa[col].to_numpy()[::1],
                label='QSSA',
                s=10,
                alpha=1,
                marker=symbol_arr[1],
                facecolors=qssa_color,
                edgecolors=qssa_color,
                linewidths=0
            )

            # Scatter plot for RCCE with hollow markers
            ax.scatter(
                x[::1],
                df_rcce[col].to_numpy()[::1],
                label='RCCE',
                s=10,
                alpha=0.9,
                marker=symbol_arr[1],
                facecolors=rcce_color,
                edgecolors=rcce_color,
                linewidths=0
            )

            # Set y-axis to log scale if necessary
            if col in range_cols:
                ax.set_yscale('log')
                ax.set_ylim(range_cols[col])
            else:
                ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
                ax.yaxis.get_major_formatter().set_powerlimits((-4, 4))

            # Set y-axis label only for the leftmost column
            if j == 0:
                ax.set_ylabel(col, fontsize=10)

            # Set x-axis limits
            ax.set_xlim(x_range)

            # Hide y-axis tick labels for all but the first column
            if j != 0:
                ax.yaxis.set_tick_params(labelleft=False)

        # Add noise column titles
        if j == 0:  # Title for the first column (Baseline)
            axs[0, j].set_title("Baseline", fontsize=10)
        else:  # Title for the other columns (Noise %)
            axs[0, j].set_title(f"{int(target_noise * 100)}%", fontsize=10)

    # Custom legend
    legend_handles = [
        Line2D([0], [0], color='black', label='VODE', linestyle='-', linewidth=1),  # VODE line
        Line2D([0], [0], marker=symbol_arr[1], color='w', label='QSSA',  # QSSA marker
                markerfacecolor=qssa_color, markeredgecolor=qssa_color, markersize=6, linestyle='None'),
        Line2D([0], [0], marker=symbol_arr[1], color='w', label='RCCE',  # RCCE marker
                markerfacecolor=rcce_color, markeredgecolor=rcce_color, markersize=6, linestyle='None')
    ]

    fig.legend(
        handles=legend_handles,
        labels=['VODE', 'QSSA', 'RCCE'],
        loc='upper center',
        bbox_to_anchor=(0.5, 0.99),
        ncol=3,
        frameon=False
    )

    # Add x-axis label only once at the bottom of the whole figure
    fig.text(0.5, 0.03, 'x (mm)', ha='center', va='center', fontsize=10)

    # Adjust spacing
    plt.subplots_adjust(wspace=0.1, hspace=0.1, top=0.89, bottom=0.1, left=0.15)

    # Save the figure
    fig_dir = "figs/NH3"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_x_overall.png"), dpi=300)
    plt.show()
    

def plot_peak_relative_error_over_noise():
    """
    Calculates and plots the peak relative error (%) of different species using QSSA and RCCE methods
    across various noise levels on the same figure. Different species are distinguished
    using different colors from color_arr, while QSSA and RCCE are differentiated by different markers and line styles.
    The y-axis is set to a logarithmic scale, and reference lines at 10%, 50%, 100%, and 1000% errors are added.
    Legends for Species are placed in the lower right corner.
    Methods are described in the plot caption using LaTeX.
    """
    # Define Method Paths
    qssa_method_path = "figs/QSSA/data/"
    rcce_method_path = "figs/RCCE/data/"
    case_name = "case_NH3_counterflow"
    
    # Define color_arr, line_arr, symbol_arr
    line_arr = ('-','--','-.',':')
    color_arr = ('k','r','b','g', 'c','m')  # Black, Red, Blue, Green, Cyan, Magenta
    symbol_arr = ('s','o','v','^','*')         # Square, Circle, Triangle Down, Triangle Up, Star
    
    # Configuration Parameters
    target_case = 'N_CF_1'
    target_noise_assemble = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.3, 0.35, 0.40, 0.45]  # Noise levels
    data_cols = ['OH', 'NO', 'NO2', 'Qdot']  # Updated Species columns
    
    # Define Paths
    base_case_name = 'case_NH3_counterflow'
    original_data_path = os.path.join("data", base_case_name, f"{target_case}.csv")
    noise_case_name = 'case_NH3_counterflow_noise'
    
    # Verify Original Data Exists
    if not os.path.exists(original_data_path):
        raise FileNotFoundError(f"Original data file not found: {original_data_path}")

    # Read Original Data
    df_original = pd.read_csv(original_data_path)
    
    # Initialize Dictionary to Store Errors
    error_dict = {col: {'QSSA': [], 'RCCE': []} for col in data_cols}
    
    # Iterate Over Each Noise Level
    for target_noise in target_noise_assemble:
        noise_target_case = f'{target_case}_{target_noise}'
        
        # Define Predicted Data Paths
        rcce_data_path = os.path.join(rcce_method_path, noise_case_name, noise_target_case, "predicted_X.csv")
        qssa_data_path = os.path.join(qssa_method_path, noise_case_name, noise_target_case, "predicted_X.csv")
        
        # Read RCCE Predicted Data
        if os.path.exists(rcce_data_path):
            df_rcce = pd.read_csv(rcce_data_path)
        else:
            print(f"Warning: RCCE data file not found: {rcce_data_path}")
            # Append NaN for RCCE if file not found
            for col in data_cols:
                error_dict[col]['RCCE'].append(float('nan'))
            df_rcce = None

        # Read QSSA Predicted Data
        if os.path.exists(qssa_data_path):
            df_qssa = pd.read_csv(qssa_data_path)
        else:
            print(f"Warning: QSSA data file not found: {qssa_data_path}")
            # Append NaN for QSSA if file not found
            for col in data_cols:
                error_dict[col]['QSSA'].append(float('nan'))
            df_qssa = None

        # Compute Peak Relative Errors if Data Exists
        for col in data_cols:
            # Initialize errors as NaN
            error_qssa = float('nan')
            error_rcce = float('nan')
            
            # Calculate Peak Relative Error for QSSA
            if df_qssa is not None and col in df_qssa.columns and col in df_original.columns:
                original_max = df_original[col].max()
                qssa_max = df_qssa[col].max()
                if original_max != 0:
                    error_qssa = np.abs((qssa_max - original_max) / original_max) * 100
                else:
                    # Handle division by zero
                    error_qssa = np.nan
                    print(f"Warning: Original max value for species '{col}' is zero. Relative error set to NaN for QSSA.")
            else:
                print(f"Warning: Missing data for QSSA method, species '{col}' at noise level {target_noise}.")
            
            # Calculate Peak Relative Error for RCCE
            if df_rcce is not None and col in df_rcce.columns and col in df_original.columns:
                original_max = df_original[col].max()
                rcce_max = df_rcce[col].max()
                if original_max != 0:
                    error_rcce = np.abs((rcce_max - original_max) / original_max) * 100
                else:
                    # Handle division by zero
                    error_rcce = np.nan
                    print(f"Warning: Original max value for species '{col}' is zero. Relative error set to NaN for RCCE.")
            else:
                print(f"Warning: Missing data for RCCE method, species '{col}' at noise level {target_noise}.")

            # Append Errors to Dictionary
            error_dict[col]['QSSA'].append(error_qssa)
            error_dict[col]['RCCE'].append(error_rcce)

    # Prepare Noise Levels for Plotting (Convert to Percentage)
    noise_levels_percent = [n * 100 for n in target_noise_assemble]

    # Initialize Plot
    plt.figure(figsize=(4.7, 4))  # Reduced figure size

    # Define Colors for Species from color_arr
    species_colors = {col: color_arr[i % len(color_arr)] for i, col in enumerate(data_cols)}
    
    # Define Markers and Line Styles for Methods
    method_styles = {
        'QSSA': {'marker': 'o', 'linestyle': '-'},    # Circle marker with solid line
        'RCCE': {'marker': 's', 'linestyle': '--'},  # Square marker with dashed line
        # Add more methods and their styles here if needed
    }

    # Plotting each species-method combination
    for col in data_cols:
        for idx, method in enumerate(['QSSA', 'RCCE']):
            errors = error_dict[col][method]
            color = species_colors.get(col, 'gray')  # Default color if not specified
            style = method_styles.get(method, {'marker': 'D', 'linestyle': '-'})  # Default styles
            marker = style['marker']
            linestyle = style['linestyle']
            
            # Assign label only once per species to avoid duplicate legend entries
            if method == 'QSSA':
                label = 'HRR' if col == 'Qdot' else col  # Change 'Qdot' to 'HRR'
            else:
                label = None  # No label for methods
            
            plt.plot(
                noise_levels_percent,
                errors,
                label=label,
                marker=marker,
                linestyle=linestyle,
                color=color,
                # markerfacecolor='none',       # Hollow markers
                markeredgecolor=color,        # Marker edge color matches line color
                markersize=6,                 # Reduced marker size
                linewidth=1,                  # Reduced line width
                alpha=0.8                     # Slight transparency for better visibility
            )

    # Set y-axis to logarithmic scale
    plt.yscale('log')

    # Add reference lines at 10%, 50%, 100%, and 1000% relative errors
    for ref_error in [10, 50, 100, 1000]:
        plt.axhline(y=ref_error, color='black', linestyle='--', linewidth=0.8)
        # Position the text slightly above the reference line and aligned to the right
        plt.text(noise_levels_percent[-1], ref_error * 1.05, f'{ref_error}%', color='black', fontsize=10,
                 verticalalignment='bottom', horizontalalignment='right')

    # Customize Plot
    plt.xlabel('Noise Level (%)', fontsize=10)
    plt.ylabel('Peak Value Reconstruction Error (%)', fontsize=10)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

    # Create Custom Legends
    # Legend for Species (Colors)
    species_handles = [
        Line2D([0], [0], marker='o', color='w', label='HRR' if species == 'Qdot' else species,
               markerfacecolor=color, markersize=10) 
        for species, color in species_colors.items()
    ]
    # Position the species legend in the lower right corner
    species_legend = plt.legend(handles=species_handles, loc='lower right', fontsize=10, title_fontsize=10)

    # Add the species legend to the plot
    plt.gca().add_artist(species_legend)

    # Adjust Layout to prevent overlap
    plt.tight_layout()

    # Create Directory for Saving Figures
    fig_dir = os.path.join("figs", "NH3", "peak_errors")
    os.makedirs(fig_dir, exist_ok=True)

    # Save the Figure
    output_path = os.path.join(fig_dir, "Max_Relative_Error_All_Species_Combined.png")
    plt.savefig(output_path, dpi=300)
    print(f"Plot saved to {output_path}")

    # Display the Plot
    plt.show()

###############MAIN##################
def main():
    plot_results_over_T()
    plot_results_over_T_overall()
    plot_results_over_space()
    # plot_results_over_space_noise()
    # plot_peak_relative_error_over_noise()


if __name__ == "__main__":
    main()
