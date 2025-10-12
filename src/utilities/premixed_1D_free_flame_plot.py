import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import pandas as pd
import numpy as np
import math
import os

def plot_temperature_against_grid(directory_path, fig_dir):
    # Create a figure for the plots
    plt.figure(figsize=(6, 4))

    # Get the list of CSV files
    csv_files = [file_name for file_name in os.listdir(directory_path) if file_name.endswith('.csv')]

    # Extract strain rates from file names and sort the files based on strain rate
    strain_rates = []
    for file_name in csv_files:
        # Extract the strain rate from the file name
        try:
            strain_rate = float(file_name.split('_')[2])
            strain_rates.append((strain_rate, file_name))
        except ValueError:
            continue  # Skip files that don't match the expected naming pattern

    # Sort the files by strain rate
    sorted_files = sorted(strain_rates, key=lambda x: x[0])

    # Create a color map
    colors = cm.jet(np.linspace(1, 0.7, len(sorted_files)))

    # Loop through all sorted files and plot them
    for idx, (strain_rate, file_name) in enumerate(sorted_files):
        file_path = os.path.join(directory_path, file_name)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        # Check if 'normalized_grid' and 'T' columns exist in the file
        if 'normalized_grid' in df.columns and 'T' in df.columns and 'CEM' in df.columns:
            # Scatter plot T against normalized_grid with a specific color
            sc = plt.scatter(df['normalized_grid'], df['T'], c=df['CEM'], cmap='jet')
    cbar = plt.colorbar(sc)
    cbar.set_label('CEM', fontsize=14)

    # Add labels and title
    plt.xlabel('x/L [-]', fontsize=14)
    plt.ylabel('Temperature [K]', fontsize=14)
    plt.legend(loc='best', fontsize='small')
    plt.grid(True)
     
    plt.xlim(0.2, 0.6)
    
    # Save the plot as a PDF file
    fig_name = 'plot_1D_counter_flow_x_T.pdf'
    fig_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_path)
    plt.show()



def plot_temperature_against_MF(directory_path, fig_dir):
    """
    Plots temperature against mixture fraction, coloring the lines based on CEM values.
    
    Parameters:
        directory_path (str): Path to the directory containing the CSV files.
    """
    # Create a figure for the plot
    plt.figure(figsize=(6, 4))

    # Get the list of CSV files
    csv_files = [file_name for file_name in os.listdir(directory_path) if file_name.endswith('.csv')]

    # Extract strain rates from file names and sort the files based on strain rate
    strain_rates = []
    for file_name in csv_files:
        # Extract the strain rate from the file name
        try:
            strain_rate = float(file_name.split('_')[2])
            strain_rates.append((strain_rate, file_name))
        except ValueError:
            continue  # Skip files that don't match the expected naming pattern

    # Sort the files by strain rate
    sorted_files = sorted(strain_rates, key=lambda x: x[0])

    # Loop through all sorted files and plot them
    for idx, (strain_rate, file_name) in enumerate(sorted_files):
        file_path = os.path.join(directory_path, file_name)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Check if 'MF', 'T', and 'CEM' columns exist in the file
        if 'MF' in df.columns and 'T' in df.columns and 'CEM' in df.columns:
            sc = plt.scatter(df['MF'], df['T'], c=df['CEM'], cmap='jet', edgecolor='none')
    
    # Add a color bar for the CEM values
    cbar = plt.colorbar(sc)
    cbar.set_label('CEM', fontsize=14)
    
    # Add labels and title
    plt.xlabel('Mixture Fraction $\\Xi$ [-]', fontsize=14)
    plt.ylabel('Temperature [K]', fontsize=14)
    plt.legend(loc='best', fontsize='small')
    plt.grid(True)
    
    # Save the plot as a PDF file
    fig_name = 'plot_1D_counter_flow_MF_T_CEM.pdf'
    fig_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_path)
    


def plot_CEM_against_temperature(directory_path, fig_dir):
    """
    Plots temperature against mixture fraction, coloring the lines based on CEM values.
    
    Parameters:
        directory_path (str): Path to the directory containing the CSV files.
    """
    # Create a figure for the plot
    plt.figure(figsize=(6, 4))

    # Get the list of CSV files
    csv_files = [file_name for file_name in os.listdir(directory_path) if file_name.endswith('.csv')]

    # Extract strain rates from file names and sort the files based on strain rate
    strain_rates = []
    for file_name in csv_files:
        # Extract the strain rate from the file name
        try:
            strain_rate = float(file_name.split('_')[2])
            strain_rates.append((strain_rate, file_name))
        except ValueError:
            continue  # Skip files that don't match the expected naming pattern

    # Sort the files by strain rate
    sorted_files = sorted(strain_rates, key=lambda x: x[0])

    # Loop through all sorted files and plot them
    for idx, (strain_rate, file_name) in enumerate(sorted_files):
        file_path = os.path.join(directory_path, file_name)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Check if 'MF', 'T', and 'CEM' columns exist in the file
        if 'CEM' in df.columns and 'T' in df.columns:
            sc = plt.scatter(df['T'], df['CEM'], c='black', edgecolor='none')
    
    # Add labels and title
    plt.xlabel('Temperature [K]', fontsize=14)
    plt.ylabel('CEM [-]', fontsize=14)
    plt.legend(loc='best', fontsize='small')
    plt.grid(True)
    
    # Save the plot as a PDF file
    fig_name = 'plot_1D_counter_flow_T_CEM.pdf'
    fig_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_path)
    
        
def plot_all_against_grid(directory_path, fig_dir, properties):
    """
    Plots multiple properties (like T, density, H2, OH, CH4) against the normalized grid
    in a square layout of subplots for all CSV files in the directory.
    
    Parameters:
        directory_path (str): Path to the directory containing the CSV files.
        properties (list): List of properties to plot (e.g., ['T', 'density', 'H2', 'OH', 'CH4']).
    """
    # Get the list of CSV files
    csv_files = [file_name for file_name in os.listdir(directory_path) if file_name.endswith('.csv')]

    # Extract strain rates from file names and sort the files based on strain rate
    strain_rates = []
    for file_name in csv_files:
        # Extract the strain rate from the file name
        try:
            strain_rate = float(file_name.split('_')[2])
            strain_rates.append((strain_rate, file_name))
        except ValueError:
            continue  # Skip files that don't match the expected naming pattern

    # Sort the files by strain rate
    sorted_files = sorted(strain_rates, key=lambda x: x[0])

    # Number of properties to plot determines the number of subplots
    n_properties = len(properties)
    
    # Calculate the number of rows and columns for the grid layout
    n_cols = math.ceil(math.sqrt(n_properties))
    n_rows = math.ceil(n_properties / n_cols)
    
    # Create a figure with subplots arranged in a grid (square) layout
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    axes = axes.flatten()  # Flatten in case it's a 2D array of axes
    
    # Create a color map
    colors = cm.jet(np.linspace(1, 0, len(sorted_files)))

    # Loop through all sorted files and plot them
    for idx, (strain_rate, file_name) in enumerate(sorted_files):
        file_path = os.path.join(directory_path, file_name)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        # Check if 'normalized_grid' exists and plot each property
        if 'normalized_grid' in df.columns:
            for i, prop in enumerate(properties):
                if prop in df.columns:
                    axes[i].plot(df['normalized_grid'], df[prop], color=colors[idx])
                    axes[i].set_ylabel(prop, fontsize=14)
                    axes[i].grid(True)

    # Set labels for the x-axis in the bottom row of subplots
    for ax in axes[-n_cols:]:  # Only set x-label for bottom row
        ax.set_xlabel('x/L', fontsize=14)
    
    # Set font size for ticks in all subplots
    for ax in axes:
        ax.tick_params(axis='both', labelsize=14)
    
    # Hide unused subplots (if there are any)
    for i in range(n_properties, len(axes)):
        fig.delaxes(axes[i])

    # Add a legend to the first subplot (optional, you can modify it to your needs)
    axes[0].legend(loc='best', fontsize='small')

    # Adjust layout and save the figure
    plt.tight_layout()
    fig_name = 'plot_all_properties.pdf'
    fig_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_path)
    plt.show()


def main():
    # Directory containing the CSV files
    # premixed 
    directory_path = 'data/case_CH4_freeflame_premixed'
    fig_dir = 'figs/free_flame'
    # non-premixed
    # directory_path = 'data/case_CH4_counterflow'
    # fig_dir = 'figs/non-premixed'
    
    # Create the fig_dir if it doesn't exist
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    
    properties_all = ['T', 'density', 'MF', 'CH4', 'CEM', 'CO2', 'CO', 'H2O', 'N2', 'O2','H2', 'OH']
    # Call the function to plot
    plot_temperature_against_grid(directory_path, fig_dir)
    plot_temperature_against_MF(directory_path, fig_dir)
    plot_CEM_against_temperature(directory_path, fig_dir)
    plot_all_against_grid(directory_path, fig_dir, properties=properties_all)

if __name__ == "__main__":
    main()

