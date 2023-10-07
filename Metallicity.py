#1 /usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification', 'Metals_IntraCluster_Stars_Mass']  # Change to your desired column names
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
    df.loc[:, 'metallicity'] = np.log10((df['Metals_IntraCluster_Stars_Mass'] / df['Intracluster_Stars_Mass']) / 0.02) #+ 9.0
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
    metals = np.log10(np.array(property_to_calculate))

    mi = 6.5
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

       w = np.where((halo >= min_bin) & (halo < max_bin) & (ihstars > 0))[0] 
      
       mean_bin = np.mean(metals[w])
       mean_values.append(mean_bin)
       
       std_bin = np.std(metals[w])
       std_values.append(std_bin)
       
       median_bin = np.median(metals[w])
       median_values.append(median_bin)

       error_mean = std_bin / np.sqrt(len(metals[w]))
       error_mean_values.append(error_mean)
       
       # print(i, min_bin, max_bin, mean_bin, std_bin) 
        
    
    return mean_values, std_values, median_values, mid_bin_values, error_mean_values

# -------------------------------------------------------------------------

def Metallicity_smass(df, property_1, property_2, property_3, titles, save_filename):
    num_plots = len(df)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 18))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    for idx, (df, title) in enumerate(zip(df, titles)):
        row = idx // num_cols
        col = idx % num_cols
        ax = axs[row, col]
        
        mass = np.log10(df[property_1] * 1.0e10 / Hubble_h)
        # print(mass)
        Z = np.log10((df[property_2] / df[property_3]) / 0.02) + 9.0
        # print(Z)

        # overplot Tremonti et al. 2003 (h=0.7)
        w = np.arange(7.0, 13.0, 0.1)
        Zobs = -1.492 + 1.847*w - 0.08026*w*w

        if 'z = 0.00' in title:
            # Conversion from Kroupa IMF to Slapeter IMF to Chabrier IMF
            ax.plot(np.log10((10**w *1.5 /1.8)), Zobs, 'b', lw=2.0, label='Tremonti et al. 2003')

        ax.scatter(mass, Z, marker='o', s=1, c='gray', alpha=0.5, label='Model galaxies')
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        ax.axis([7.5, 12.0, 7.0, 9.5])
        
        leg = ax.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

     # Blank out remaining unused subplots
    for idx in range(num_plots, num_cols * num_rows):
        row = idx // num_cols
        col = idx % num_cols
        axs[row, col].axis('off')

    plt.tight_layout()
    save_plot(save_filename)

def Metallicity_hmass(df, property_1, property_2, property_3, titles, save_filename):
    num_plots = len(df)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 18))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    for idx, (df, title) in enumerate(zip(df, titles)):
        row = idx // num_cols
        col = idx % num_cols
        ax = axs[row, col]
        
        mass = np.log10(df[property_1] * 1.0e10 / Hubble_h)
        # print(mass)
        Z = np.log10((df[property_2] / df[property_3]) / 0.02) #+ 9.0
        # print(Z)

        ax.scatter(mass, Z, marker='o', s=1, c='gray', alpha=0.5, label='Model galaxies')
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        ax.axis([11.0, 15.5, 5.0, 9.5])
        
        leg = ax.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

     # Blank out remaining unused subplots
    for idx in range(num_plots, num_cols * num_rows):
        row = idx // num_cols
        col = idx % num_cols
        axs[row, col].axis('off')

    plt.tight_layout()
    save_plot(save_filename)          
    

def Metallicity_ihs(df, property_1, property_2, property_3, titles, save_filename):
    num_plots = len(df)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 18))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    for idx, (df, title) in enumerate(zip(df, titles)):
        row = idx // num_cols
        col = idx % num_cols
        ax = axs[row, col]
        
        mass = np.log10(df[property_1] * 1.0e10 / Hubble_h)
        # print(mass)
        Z = np.log10((df[property_2] / df[property_3]) / 0.02) + 9.0
        # print(Z)

        ax.scatter(mass, Z, marker='o', s=1, c='gray', alpha=0.5, label='Model galaxies')
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        ax.axis([7.5, 12.0, 5.0, 9.5])
        
        leg = ax.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

     # Blank out remaining unused subplots
    for idx in range(num_plots, num_cols * num_rows):
        row = idx // num_cols
        col = idx % num_cols
        axs[row, col].axis('off')

    plt.tight_layout()
    save_plot(save_filename)  

def Metallicity_ihsfraction(df, property_1, property_2, property_3, titles, save_filename):
    num_plots = len(df)
    num_cols = 2
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 18))
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    
    for idx, (df, title) in enumerate(zip(df, titles)):
        row = idx // num_cols
        col = idx % num_cols
        ax = axs[row, col]
        
        mass = df[property_1]
        # print(mass)
        Z = np.log10((df[property_2] / df[property_3]) / 0.02) + 9.0
        # print(Z)

        ax.scatter(mass, Z, marker='o', s=1, c='gray', alpha=0.5, label='Model galaxies')
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
        ax.set_xlabel(r'IHS Fraction')  # and the x-axis labels
            
        # Set the x and y axis minor ticks
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
            
        ax.axis([0, 0.15, 5.0, 9.5])
        
        leg = ax.legend(loc='lower right')
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

     # Blank out remaining unused subplots
    for idx in range(num_plots, num_cols * num_rows):
        row = idx // num_cols
        col = idx % num_cols
        axs[row, col].axis('off')

    plt.tight_layout()
    save_plot(save_filename)       

def Metallicity_all(save_filename):

    csv_files = '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4455.0.csv'
    titles = 'z = 0.00'

    print('Processing simulation data...', titles)
    df = pd.read_csv(csv_files)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_diluted = dilute_dataframe(df_filtered, target_rows)
    df_calculated = perform_calculations(df_diluted)


    mean_values, std_values, median_values, mid_bin_values, error_mean_values = calculate_statistics(df, df_calculated['hmass'], df_calculated['metallicity'])
        
    mean_values = np.array(mean_values) 
    std_values = np.array(std_values)
    median_values = np.array(median_values)
    error_mean_values = np.array(error_mean_values)
        
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    
    axs[0, 0].scatter(np.log10(df_calculated['hmass']),df_calculated['metallicity'],s=1,c='gray',alpha=0.5)
    # axs[0, 0].set_title('Plot 1: sin(x)')
    axs[0, 0].set_ylim(-3,1)
    axs[0, 0].plot(mid_bin_values, mean_values, color='blue', label='Mean')
    axs[0, 0].plot(mid_bin_values, median_values, color='green', label='Median', linestyle = '--')
    axs[0, 0].fill_between(mid_bin_values, (mean_values-std_values), (mean_values+std_values), color='orange', alpha=0.3, label=r'1$\sigma$ Std. Dev.')
    axs[0, 0].fill_between(mid_bin_values, (mean_values-error_mean_values), (mean_values+error_mean_values), color='red', alpha=0.3, label=r'Error in mean')
    axs[0, 0].set_ylabel(r'$\log_{10}[Z_{IHS}/Z_{\odot}]$')  # Set the y...
    axs[0, 0].set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[0, 1].scatter(np.log10(df_calculated['smass']),df_calculated['metallicity'],s=1,c='gray',alpha=0.5)
    # axs[0, 1].set_title('Plot 2: cos(x)')
    w = np.arange(7.0, 13.0, 0.1)
    Zobs = -1.492 + 1.847*w - 0.08026*w*w
    axs[0, 1].set_ylim(-3,1)
    # axs[0, 1].plot(np.log10((10**w *1.5 /1.8)), Zobs, 'b', lw=2.0, label='Tremonti et al. 2003')
    axs[0, 1].set_ylabel(r'$\log_{10}[Z_{IHS}/Z_{\odot}]$')  # Set the y...
    axs[0, 1].set_xlabel(r'$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[1, 0].scatter(np.log10(df_calculated['IHM']),df_calculated['metallicity'],s=1,c='gray',alpha=0.5)
    # axs[1, 0].set_title('Plot 3: tan(x)')
    axs[1, 0].set_ylim(-3,1)
    axs[1, 0].set_ylabel(r'$\log_{10}[Z_{IHS}/Z_{\odot}]$')  # Set the y...
    axs[1, 0].set_xlabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[1, 1].scatter(df_calculated['IHM_Fraction'],df_calculated['metallicity'],s=1,c='gray',alpha=0.5)
    # axs[1, 1].set_title('Plot 4: exp(x)')
    axs[1, 1].set_ylim(-3,1)
    axs[1, 1].set_ylabel(r'$\log_{10}[Z_{IHS}/Z_{\odot}]$')  # Set the y...
    axs[1, 1].set_xlabel('Intrahalo Stars Fraction')  # and the x-axis labels
    
    plt.tight_layout()
    save_plot(save_filename)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------    
# -------------------------------------------------------------------------

# # List of CSV files to process along with corresponding titles
# csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4455.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4456.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4457.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4459.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4460.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4461.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4462.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4463.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4464.0.csv',
#              '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4465.0.csv']  # Add your CSV filenames here
# titles = ['z = 0.00', 'z = 0.508', 'z = 1.0', 'z = 2.07', 'z = 3.06', 'z = 4.17', 'z = 5.288', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles
# # Initialize lists to store DataFrames for old and new datasets
# datasets = []

# # Loop through the specified CSV files for 1 dataset
# for idx, filename in enumerate(csv_files):
#     print('Processing simulation data...', titles[idx])
#     df = pd.read_csv(filename)
#     df = df[columns_to_extract]
#     df_filtered = df[df['Galaxy_Classification'] == 0]
#     df_diluted = dilute_dataframe(df_filtered, target_rows)
#     df_calculated = perform_calculations(df_diluted)
#     datasets.append(df_calculated)

# Metallicity_smass(datasets, 'Total_Stellar_Mass', 'Metals_IntraCluster_Stars_Mass', 'Intracluster_Stars_Mass', titles, 'Metallicity_smass.png')
# Metallicity_hmass(datasets, 'Mvir', 'Metals_IntraCluster_Stars_Mass', 'Intracluster_Stars_Mass', titles, 'Metallicity_hmass.png')
# Metallicity_ihs(datasets, 'Intracluster_Stars_Mass', 'Metals_IntraCluster_Stars_Mass', 'Intracluster_Stars_Mass', titles, 'Metallicity_ihs.png')
# Metallicity_ihsfraction(datasets, 'IHM_Fraction', 'Metals_IntraCluster_Stars_Mass', 'Intracluster_Stars_Mass', titles, 'Metallicity_ihsfraction.png')
Metallicity_all('Metallicity_all.png')