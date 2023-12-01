#1 /usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification']  # Change to your desired column names
target_rows = 200000
Hubble_h = 0.73
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
    df.loc[:, 'halo_mass'] = df['Mvir'] * 1.0e10 / Hubble_h
    df.loc[:, 'IHM_Fraction'] = df['Intracluster_Stars_Mass'] / (0.17 * df['Mvir'])  # Modify this calculation as needed
    return df

# Function to save the plot with higher quality
def save_plot(save_filename):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, '/Users/michaelbradley/Documents/Honours/Images/', save_filename)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.show()

# -------------------------------------------------------------------------
    
def plot_IHMFraction_vs_redshift(df_list, property_1, property_2, titles, save_filename):
    plt.figure()
    
    colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
    redshifts = [0.00, 0.0199, 0.0414, 0.0645, 0.0893, 0.1159, 0.1749, 0.2075, 0.2798]
    haloes = []
    ihsmf = []

    for df, color, title in zip(df_list, colors, titles):
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))

        w = np.where((halo >= 14) & (ihs_frac > 0))[0]
        hmass = halo[w]
        ihs = ihs_frac[w]
        ihs_perc = ihs * 100

        haloes.append(hmass)
        ihsmf.append(ihs_perc)

    fraction = ihsmf
    hmass_ = haloes

    plt.scatter(hmass_[0],fraction[0], alpha=0.7, c='white', edgecolors='black', s = 70, marker = 'h')
    # plt.scatter(hmass_[1],fraction[1], s=1, alpha=0.6, c='k', label='z = 0.0199')
    # plt.scatter(hmass_[2],fraction[2], s=1, alpha=0.6, c='k', label='z = 0.0414')
    # plt.scatter(hmass_[3],fraction[3], s=1, alpha=0.6, c='k', label='z = 0.0645')
    # plt.scatter(hmass_[4],fraction[4], s=1, alpha=0.6, c='k', label='z = 0.0893')
    # plt.scatter(hmass_[4],fraction[4], s=1, alpha=0.4, c='k', label='z = 0.1159')
    
    plt.xlim(14,15.2)
    plt.ylim(0, 20)
    plt.xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')
    plt.ylabel('Intrahalo stellar mass fraction (%)')

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
# -------------------------------------------------------------------------

# List of CSV files to process along with corresponding titles
csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4406.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4427.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4407.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4428.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4408.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4409.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4410.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4411.0.csv',]  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 0.1159',
           'z = 0.1749', 'z = 0.2075',
            'z = 0.2798']  # Add corresponding titles

datasets = []

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    df_diluted = dilute_dataframe(df_calculated, target_rows)
    datasets.append(df_diluted)

print('Intrahalo mass fraction vs. Redshift')
plot_IHMFraction_vs_redshift(datasets, 'halo_mass', 'IHM_Fraction', titles, 'IHM_Fraction_zoom_scatter.png')
