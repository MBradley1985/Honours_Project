import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification']  # Change to your desired column names
target_rows = 7500

# -------------------------------------------------------------------------

# Function to dilute a dataset to the target number of rows
def dilute_dataframe(df, target_rows):
    df = df[df['Galaxy_Classification'] == 0]
    if len(df) > target_rows:
        return df.sample(n=target_rows, random_state=42)
    else:
        return df

# Function to save the plot with higher quality
def save_plot(save_filename):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, '/Users/michaelbradley/Documents/Honours/Images/', save_filename)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.show()

# -------------------------------------------------------------------------

def divide_and_store_custom_mass_ranges(df, bin_edges):
    adjusted_mvir = np.log10(df['Mvir'] * 1.0e10 / 0.73)
    mass_ranges = pd.cut(adjusted_mvir, bins=bin_edges)
    grouped = df.groupby(mass_ranges)['Intracluster_Stars_Mass'].apply(list)
    return grouped

# -------------------------------------------------------------------------

def process_subsets(data_per_file, datasets, property_1, subset_titles, main_titles, rows, cols, save_filename):
    BoxSize = 500.0  # millennium
    box_fraction = 1.0
    
    volume = (BoxSize / 0.73)**3.0 * box_fraction
    
    mi = 7.5
    ma = 12.75
    binwidth = 0.25
    NB = int((ma - mi) / binwidth)
    
    num_files = len(data_per_file)
    subset_colors = ['red', 'green', 'blue']  # List of colors for subsets
    legend_labels = [r'$10^{10.5} < M_{\mathrm{halo}} < 10^{12}$',
                 r'$10^{12} < M_{\mathrm{halo}} < 10^{13.5}$',
                 r'$10^{13.5} < M_{\mathrm{halo}}$']  # List of custom legend labels for the three bins

    
    fig, axes = plt.subplots(rows, cols, figsize=(12, 18))
    
    custom_lines = []  # For creating a custom legend
    
    for file_idx, (subset_data, df, subset_title, main_title) in enumerate(zip(data_per_file, datasets, subset_titles, main_titles)):
        row = file_idx // cols
        col = file_idx % cols
        ax = axes[row, col]
        
        for data, color, label in zip(subset_data, subset_colors, legend_labels):
            (counts, binedges) = np.histogram(np.log10(data * 1.0e10 / 0.73), range=(mi, ma), bins=NB)
            xaxeshisto = binedges[:-1] + 0.5 * binwidth
            line, = ax.plot(xaxeshisto, counts / volume / binwidth, color=color, label=label, alpha=0.5)
            ax.fill_between(xaxeshisto, 0, counts / volume / binwidth, color=color, alpha=0.1)
            custom_lines.append(line)
        
        (counts, binedges) = np.histogram(np.log10(df[property_1] * 1.0e10 / 0.73), range=(mi, ma), bins=NB)
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        line, = ax.plot(xaxeshisto, counts / volume / binwidth, color='black', label='Overall')
        custom_lines.append(line)
        
        ax.set_yscale('log')
        ax.set_ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$')
        ax.set_title(main_title)
        
        handles_1 = custom_lines[:3]  # Assuming the first three lines correspond to your custom legend labels
        
        handles_2 = [Line2D([], [], color='black', label='Overall')]  # Line for the second histogram dataset
        
        ax.legend(handles=handles_1 + handles_2, loc='upper right', frameon=False, fontsize='small')

    # Hide unused subplots
    for file_idx in range(num_files, rows * cols):
        row = file_idx // cols
        col = file_idx % cols
        axes[row, col].axis('off')
    
    plt.tight_layout()
    save_plot(save_filename)   

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4406.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4427.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4407.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4428.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4420.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4422.0.csv']  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 1.0', 'z = 2.07']  # Add corresponding titles
titles_sub = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 1.0', 'z = 2.07']  # Add corresponding titles

# List of CSV files to process along with corresponding titles
csv_files2 = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4410.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4414.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4419.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4420.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4421.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4422.0.csv']  # Add your CSV filenames here
titles2 = ['z = 0.00', 'z = 0.201', 'z = 0.508', 'z = 0.827', 'z = 1.077', 'z = 1.503', 'z = 2.07']  # Add corresponding titles
titles2_sub = ['z = 0.00', 'z = 0.201', 'z = 0.508', 'z = 0.827', 'z = 1.077', 'z = 1.503', 'z = 2.07']  # Add corresponding titles

# List of CSV files to process along with corresponding titles
csv_files3 = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4420.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4422.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4423.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4424.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4425.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4426.0.csv']  # Add your CSV filenames here
titles3 = ['z = 0.00', 'z = 1.0', 'z = 2.07', 'z = 4.17', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles
titles3_sub = ['z = 0.00', 'z = 1.0', 'z = 2.07', 'z = 4.17', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles

datasets = []
datasets2 = []
datasets3 = []
grouped_data_per_file = []
grouped_data_per_file2 = []
grouped_data_per_file3 = []
# datasets_list = [grouped_data_per_file,datasets3]  # List of datasets

custom_bin_edges = [10.5, 12, 13.5, 17]  # Modify these edges as needed

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_diluted = dilute_dataframe(df, target_rows)
    datasets.append(df_diluted)
    mass_range_data = divide_and_store_custom_mass_ranges(df_diluted, custom_bin_edges)
    mass_range_arrays = [np.array(data) for data in mass_range_data]
    grouped_data_per_file.append(mass_range_arrays)
    
for idx, filename in enumerate(csv_files2):
    print('Processing simulation data...', titles2[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_diluted = dilute_dataframe(df, target_rows)
    datasets2.append(df_diluted)
    mass_range_data = divide_and_store_custom_mass_ranges(df_diluted, custom_bin_edges)
    mass_range_arrays = [np.array(data) for data in mass_range_data]
    grouped_data_per_file2.append(mass_range_arrays)

for idx, filename in enumerate(csv_files3):
    print('Processing simulation data...', titles3[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_diluted = dilute_dataframe(df, target_rows)
    datasets3.append(df_diluted)
    mass_range_data = divide_and_store_custom_mass_ranges(df_diluted, custom_bin_edges)
    mass_range_arrays = [np.array(data) for data in mass_range_data]
    grouped_data_per_file3.append(mass_range_arrays)

process_subsets(grouped_data_per_file, datasets, columns_to_extract[1], titles_sub, titles, 4, 2, 'IHMfunction_1.png')
process_subsets(grouped_data_per_file2, datasets2, columns_to_extract[1], titles2_sub, titles2, 4, 2, 'IHMfunction_2.png')
process_subsets(grouped_data_per_file3, datasets3, columns_to_extract[1], titles3_sub, titles3, 4, 2, 'IHMfunction_3.png')