import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import matplotlib.ticker as ticker

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification']  # Change to your desired column names
target_rows = 10000
Hubble_h = 0.73
binwidth = 0.25

# -------------------------------------------------------------------------

# Function to dilute a dataset to the target number of rows
def dilute_dataframe(df, target_rows):
    
    if len(df) > target_rows:
        return df.sample(n=target_rows, random_state=42)
    else:
        return df

# Function to perform calculations on columns and return a modified DataFrame
def perform_calculations(df):
    
    df = df.copy ()
    df.loc[:, 'IHM_Fraction'] = df['Intracluster_Stars_Mass'] / (0.17 * df['Mvir'])  # Modify this calculation as needed
    df.loc[:, 'IHM'] = df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h
    df.loc[:, 'hmass'] = df['Mvir'] * 1.0e10 / Hubble_h
    df.loc[:, 'smass'] = df['Total_Stellar_Mass'] * 1.0e10 / Hubble_h
    return df

# Function to save the plot with higher quality
def save_plot(save_filename):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, '/Users/michaelbradley/Documents/Honours/Images/', save_filename)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.show()

def calculate_statistics(df, property_to_bin, property_to_calculate):
    
    halo = np.log10(np.array(property_to_bin))
    ihsfraction = np.array(property_to_calculate)

    mi = 11.0
    ma = 16.5
    bin_edges = np.arange(mi, ma + binwidth, binwidth)
                   
    mean_values = []
    std_values = []
    median_values = []
    mid_bin_values = []
    error_mean_values = []

    for i in range(len(bin_edges)):

       min_bin = mi + binwidth*i
       max_bin = mi + binwidth*(i+1)
       mid_bin_values.append( mi + binwidth*(i+0.5) )

       w = np.where((halo >= min_bin) & (halo < max_bin) & (ihsfraction > 0))[0] 
      
       mean_bin = np.mean(ihsfraction[w])
       mean_values.append(mean_bin)
       
       std_bin = np.std(ihsfraction[w])
       std_values.append(std_bin)
       
       median_bin = np.median(ihsfraction[w])
       median_values.append(median_bin)

       error_mean = std_bin / np.sqrt(len(ihsfraction[w]))
       error_mean_values.append(error_mean)
       
       # print(i, min_bin, max_bin, mean_bin, std_bin) 
        
    
    return mean_values, std_values, median_values, mid_bin_values, error_mean_values

# -------------------------------------------------------------------------

def Scatter_Plot(df, property_1, property_2, titles, ax):
    
    
    # dilute_dataframe(df, target_rows)
    
    for df, title in zip(df, titles):
        df_dilute = dilute_dataframe(df, target_rows)
        ax.scatter(np.log10(df_dilute[property_1]), (df_dilute[property_2]), alpha=0.6, s=0.8, color= 'gray')

    
# -------------------------------------------------------------------------

def plot_IHMFraction_vs_hmass(df, property_1, property_2, titles, save_filename):
    num_plots = len(df)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 18))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    for idx, (df, title) in enumerate(zip(df, titles)):
        row = idx // num_cols
        col = idx % num_cols
        ax = axs[row, col]
        
        Scatter_Plot([df], property_1, property_2, [title], ax)
        
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')
        ax.set_ylabel(r'$F_{\mathrm{IHS}}\ (M_{\odot}) \left(\frac{M_{\mathrm{IHS}}}{M_{\mathrm{halo}}} \times f_{\mathrm{cosmic\ baryons}}\right)$')
        # ax.set_yscale('log')
        
        # Manually adding the title in the upper left corner
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')

        mean_values, std_values, median_values, mid_bin_values, error_mean_values = calculate_statistics(df, df['hmass'], df['IHM_Fraction'])
        
        mean_values = np.array(mean_values) 
        std_values = np.array(std_values)
        median_values = np.array(median_values)
        error_mean_values = np.array(error_mean_values)
            
        ax.xaxis.set_minor_locator(ticker.FixedLocator(mid_bin_values))
        
        ax.plot(mid_bin_values, mean_values, color='blue', label='Mean')
        ax.plot(mid_bin_values, median_values, color='green', label='Median', linestyle = '--')
        ax.fill_between(mid_bin_values, (mean_values-std_values), (mean_values+std_values), color='orange', alpha=0.3, label=r'1$\sigma$ Std. Dev.')
        ax.fill_between(mid_bin_values, (mean_values-error_mean_values), (mean_values+error_mean_values), color='red', alpha=0.3, label=r'Error in mean')
        
        # Create a custom legend with larger symbols
        leg = ax.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('medium')
      
        
     # Blank out remaining unused subplots
    for idx in range(num_plots, num_cols * num_rows):
        row = idx // num_cols
        col = idx % num_cols
        axs[row, col].axis('off')
    
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def IHM_Fraction(save_filename):

    csv_files = '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv'
    titles = 'z = 0.00'

    print('Processing simulation data...', titles)
    df = pd.read_csv(csv_files)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    df_diluted = dilute_dataframe(df_calculated, target_rows)

    mean_values, std_values, median_values, mid_bin_values, error_mean_values = calculate_statistics(df, df_calculated['hmass'], df_calculated['IHM_Fraction'])
        
    mean_values = np.array(mean_values) 
    std_values = np.array(std_values)
    median_values = np.array(median_values)
    error_mean_values = np.array(error_mean_values)

    plt.figure()

    plt.scatter(np.log10(df_diluted['hmass']), df_diluted['IHM_Fraction'], alpha=0.1, s=0.8, color='gray')
    plt.plot(mid_bin_values, mean_values, color='blue', label='Mean')
    plt.plot(mid_bin_values, median_values, color='green', label='Median', linestyle = '--')
    plt.fill_between(mid_bin_values, (mean_values-std_values), (mean_values+std_values), color='orange', alpha=0.3, label=r'1$\sigma$ Std. Dev.')
    plt.fill_between(mid_bin_values, (mean_values-error_mean_values), (mean_values+error_mean_values), color='red', alpha=0.3, label=r'Error in mean')
    
    plt.xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')
    plt.ylabel(r'$F_{\mathrm{IHS}}\ (M_{\odot}) \left(\frac{M_{\mathrm{IHS}}}{M_{\mathrm{halo}}} \times f_{\mathrm{cosmic\ baryons}}\right)$')
        
    # Create a custom legend with larger symbols
    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')
      
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------    
# -------------------------------------------------------------------------

# List of CSV files to process along with corresponding titles
csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
            '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4456.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4420.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4422.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4423.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4424.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4425.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4426.0.csv']  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.508', 'z = 1.0', 'z = 2.07', 'z = 4.17', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles


# Initialize lists to store DataFrames for old and new datasets
datasets = []

# Loop through the specified CSV files for 1 dataset
for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    datasets.append(df_calculated)
    
print('Intrahalo mass fraction vs. Halo mass')
plot_IHMFraction_vs_hmass(datasets, 'hmass', 'IHM_Fraction', titles, 'IHMFraction_vs_hmass_grid.png')

IHM_Fraction('IHM_Fraction_hmass_z0.png')