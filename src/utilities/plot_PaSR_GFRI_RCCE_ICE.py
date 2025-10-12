import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

inch = 2.54
x_range = [0, 0.2]
case_name = "PaSR"
line_arr = ('-','--','-.',':')
color_arr = ('k','r','b','y','g', 'c','m')
symbol_arr = ('s','o','v','^','*')




###########Figure 0: Plot the original data ##############
def plot_T_over_MF(target_case='KerM03_NH3CEM'):
    original_data_path = os.path.join("data", case_name, target_case) + ".csv"
    df_original = pd.read_csv(original_data_path)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(6, 6))  # Fix figure size and axes
    ax.scatter(df_original['MF'].to_numpy(), df_original['T'].to_numpy(),  # Convert to NumPy array for T
               label=f'{target_case}', 
               s=10, 
               alpha=0.8, 
               marker='o', 
               color='black')  # Original data in black
    
    # Add labels and title
    ax.set_xlabel('Mass Fraction (MF)', fontsize=12)
    ax.set_ylabel('Temperature (T)', fontsize=12)
    ax.set_title(f'Temperature vs Mass Fraction for {target_case}', fontsize=14)
    
    # Add grid and legend
    ax.grid(True)
    ax.legend()
    
    # Tight layout for better spacing
    plt.tight_layout()

    # Save the figure
    fig_dir = "figs/NH3_PaSR"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, f"1D_CEQ_{target_case}_T_vs_MF.png"), dpi=300)
    
    # Show the plot
    plt.show()


###########Figure 1: compare the QSSA and ILDM over T overall##############
# Define target cases
def plot_results_over_T(target_case, data_cols):
    
    qssa_method_path = "figs/QSSA/data/"
    rcce_method_path = "figs/RCCE/data/"
    
    
    # Initialize subplots
    fig, axs = plt.subplots(3, 3, figsize=(8, 7), sharex=True)

    # Define color for QSSA and RCCE scatter points
    qssa_color = 'blue'
    rcce_color = 'red'

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
                label=f'{target_case} - Original', 
                s=10, 
                alpha=0.3, 
                marker='o', 
                color='black')  # Original data in black
        
        # Scatter plot for QSSA
        ax.scatter(T, 
                df_qssa[col].to_numpy(),  # Convert to NumPy array
                label=f'QSSA - {target_case}', 
                s=15, 
                alpha=0.3, 
                marker='o', 
                color=qssa_color)  # QSSA data in blue (scatter points)

        # Scatter plot for RCCE
        ax.scatter(T, 
                df_rcce[col].to_numpy(),  # Convert to NumPy array
                label=f'RCCE - {target_case}', 
                s=15, 
                alpha=0.3, 
                marker='o', 
                color=rcce_color)  # RCCE data in red (scatter points)

        # Log scale for Qdot (HRR) only
        if col == 'Qdot':
            ax.set_yscale('log')
            ax.set_ylabel('HRR (J/kg/s)', fontsize=10)  # Update label for Qdot

        # Formatting
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        ax.yaxis.get_major_formatter().set_powerlimits((-2, 2))
        ax.set_ylabel(col, fontsize=10)
        # ax.set_xlim([1000, 2500])
        
        if row == 2:
            ax.set_xlabel('T(K)', fontsize=10)

    # Remove extra subplots if there are fewer data columns
    for j in range(len(data_cols), 9):
        fig.delaxes(axs[j // 3, j % 3])

    # Tight layout for better spacing
    plt.tight_layout()

    # Save the figure
    fig_dir = "figs/NH3_PaSR"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "1D_CEQ_NH3_combine_Temperature.png"), dpi=300)
    plt.show()


def plot_regression_Tcolor(target_case, data_cols):
    case_name = "PaSR"
    qssa_method_path = "figs/QSSA/data/"
    rcce_method_path = "figs/RCCE/data/"

    # Initialize subplots for 2 rows (QSSA, RCCE) and len(data_cols) columns
    fig, axs = plt.subplots(2, len(data_cols), figsize=(3*len(data_cols), 6))

    original_data_path = os.path.join("data", case_name, target_case) + ".csv"
    rcce_data_path = os.path.join(rcce_method_path, case_name, target_case, "predicted_X.csv")
    qssa_data_path = os.path.join(qssa_method_path, case_name, target_case, "predicted_X.csv")

    # Read data
    df_original = pd.read_csv(original_data_path)
    df_rcce = pd.read_csv(rcce_data_path)
    df_qssa = pd.read_csv(qssa_data_path)

    # Ensure MF and T are NumPy arrays
    MF = df_original['MF'].to_numpy()
    T = df_original['T'].to_numpy()  # Extract temperature for coloring

    # Plot each variable in its corresponding subplot
    for i, col in enumerate(data_cols):
        # QSSA subplot (row 0)
        ax = axs[0, i]
        sc = ax.scatter(df_original[col], 
                        df_qssa[col],  # Compare QSSA against original
                        c=T,           # Color by temperature
                        cmap='viridis',  # Color map for temperature
                        label=f'QSSA - {target_case}', 
                        s=15, 
                        alpha=0.8, 
                        marker='x')  # Removed the `color=qssa_color` argument
        ax.set_ylabel(f'QSSA - {col}', fontsize=10)
        
        # Calculate the range for y=x line based on the data
        x_min, x_max = df_original[col].min(), df_original[col].max()

        # Set the range for the y=x line to cover the min/max of both x and y
        min_range = x_min
        max_range = x_max

        # Plot the y=x line (diagonal) with dynamic range
        ax.plot([min_range, max_range], [min_range, max_range], color='gray', linestyle='--', label='y=x')

        # Set the x and y axis limits to be the same based on df_original's range
        ax.set_xlim([min_range, max_range])
        ax.set_ylim([min_range, max_range])

        # RCCE subplot (row 1)
        ax = axs[1, i]
        sc = ax.scatter(df_original[col], 
                        df_rcce[col],  # Compare RCCE against original
                        c=T,           # Color by temperature
                        cmap='viridis',  # Color map for temperature
                        label=f'RCCE - {target_case}', 
                        s=15, 
                        alpha=0.8, 
                        marker='x')  # Removed the `color=rcce_color` argument
        ax.set_ylabel(f'RCCE - {col}', fontsize=10)

        # Calculate the range for y=x line based on the data
        x_min, x_max = df_original[col].min(), df_original[col].max()

        # Set the range for the y=x line to cover the min/max of both x and y
        min_range = x_min
        max_range = x_max

        # Plot the y=x line (diagonal) with dynamic range
        ax.plot([min_range, max_range], [min_range, max_range], color='gray', linestyle='--', label='y=x')

        # Set the x and y axis limits to be the same based on df_original's range
        ax.set_xlim([min_range, max_range])
        ax.set_ylim([min_range, max_range])

        # Formatting for both x and y axes
        for row in range(2):  # Iterate through both rows
            axs[row, i].xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
            axs[row, i].yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
            # axs[row, i].yaxis.get_major_formatter().set_powerlimits((-3, 3))

    # Tight layout for better spacing
    plt.tight_layout()

    # Save the figure
    fig_dir = "figs/NH3_PaSR"
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(os.path.join(fig_dir, "QSSA_RCCE_vs_Original_MF.png"), dpi=300)
    plt.show()


###############MAIN##################
def main():
    target_case = 'KerM01_NH3_1e-1CEM'
    data_cols = ['OH', 'NO', 'Qdot','CEM','MF']
    plot_T_over_MF(target_case)
    plot_results_over_T(target_case, data_cols)
    plot_regression_Tcolor(target_case, data_cols) 
    # data_cols = ['NH2', 'NH', 'HNO', 'HONO', 'H2O2', 'H2NN', 'N2H2', 'NO3', 'Qdot']


if __name__ == "__main__":
    main()