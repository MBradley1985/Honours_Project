#1 /usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification']  # Change to your desired column names
target_rows = 10000
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
    ihs_non_percent1 = []
    ihs_non_percent2 = []

    for df, color, title in zip(df_list, colors, titles):
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))
        w = np.where((halo >= 14) & (ihs_frac > 0))[0]
        hmass = halo[w]
        ihs = ihs_frac[w]
        ihs_perc = ihs * 100
        w1 = np.where(ihs_perc <= 2)[0]
        w2 = np.where(ihs_perc <= 1)[0]
        ihs_non_1 = ihs_perc[w2]
        ihs_non_2 = ihs_perc[w1]
        ihs_non_perc_1 = ((len(ihs_non_1)) / (len(hmass))) * 100
        ihs_non_perc_2 = ((len(ihs_non_2)) / (len(hmass))) * 100
   
        ihs_non_percent1.append(ihs_non_perc_1)
        ihs_non_percent2.append(ihs_non_perc_2)

    # ihs_non = ihs_non_values
    ihs_non_perc1 = ihs_non_percent1
    ihs_non_perc2 = ihs_non_percent2

    plt.plot(redshifts, ihs_non_perc1, label = 'Intrahalo stars fraction < 1%',  zorder = 10)
    plt.plot(redshifts, ihs_non_perc2, label = 'Intrahalo stars fraction < 2%', zorder = 9)
    
    plt.xlim(0, 0.28)
    plt.xlabel('Redshift')
    plt.ylabel('Percent with a negligible IHM fraction')

    # Create a custom legend with larger symbols
    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')
      
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def plot_IHMFraction_vs_redshift_2(df_list, property_1, property_2, titles, save_filename):
    plt.figure()
    
    colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
    redshifts = [0.00, 0.0199, 0.0414, 0.0645, 0.0893, 0.1159, 0.1749, 0.2075, 0.2798]
    ihs_non_percent1 = []
    ihs_non_percent2 = []
    ihs_non_percent3 = []


    for df, color, title in zip(df_list, colors, titles):
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))


        a = np.where((halo >= 10.5) & (halo < 12) & (ihs_frac > 0))[0]
        hmass = halo[a]
        ihs = ihs_frac[a]
        ihs_perc = ihs * 100

        b = np.where((halo >= 12) & (halo < 13.5)  & (ihs_frac > 0))[0]
        hmass2 = halo[b]
        ihs2 = ihs_frac[b]
        ihs_perc2 = ihs2 * 100

        c = np.where((halo >= 13.5) & (ihs_frac > 0))[0]
        hmass3 = halo[c]
        ihs3 = ihs_frac[c]
        ihs_perc3 = ihs3 * 100

        w = np.where(ihs_perc >= 4)[0]
        w2 = np.where(ihs_perc2 >= 4)[0]
        w3 = np.where(ihs_perc3 >= 4)[0]
        
        ihs_non_1 = ihs_perc[w]
        ihs_non_perc_1 = ((len(ihs_non_1)) / (len(hmass))) * 100
        ihs_non_percent1.append(ihs_non_perc_1)
        print(len(hmass))
        print(len(ihs_non_1))
        print(ihs_non_perc_1)

        ihs_non_2 = ihs_perc2[w2]
        ihs_non_perc_2 = ((len(ihs_non_2)) / (len(hmass2))) * 100
        ihs_non_percent2.append(ihs_non_perc_2)

        ihs_non_3 = ihs_perc3[w3]
        ihs_non_perc_3 = ((len(ihs_non_3)) / (len(hmass3))) * 100
        ihs_non_percent3.append(ihs_non_perc_3)

    # ihs_non = ihs_non_values
    ihs_non_perc1 = ihs_non_percent1
    ihs_non_perc2 = ihs_non_percent2
    ihs_non_perc3 = ihs_non_percent3


    plt.plot(redshifts, ihs_non_perc1, label = r'$10^{10.5} < M_{\mathrm{halo}} < 10^{12}$', c = 'r', zorder = 10)
    plt.plot(redshifts, ihs_non_perc2, label = r'$10^{12} < M_{\mathrm{halo}} < 10^{13.5}$', c = 'g', zorder = 9)
    plt.plot(redshifts, ihs_non_perc3, label = r'$10^{13.5} < M_{\mathrm{halo}}$', c = 'b', zorder = 8)

    
    plt.xlim(0, 0.28)
    plt.xlabel('Redshift')
    plt.ylabel('Percent with measurable IHM fractions')

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
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4406.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4427.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4407.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4428.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4408.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4409.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4410.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4411.0.csv',]  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 0.1159',
           'z = 0.1749', 'z = 0.2075', 'z = 0.2798']  # Add corresponding titles

datasets = []

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    datasets.append(df_calculated)

# print('Intrahalo mass fraction vs. Redshift')
# plot_IHMFraction_vs_redshift(datasets, 'halo_mass', 'IHM_Fraction', titles, 'Neg_IHM_Fraction_zoom.png')

print('Intrahalo mass fraction vs. Redshift')
plot_IHMFraction_vs_redshift_2(datasets, 'halo_mass', 'IHM_Fraction', titles, 'Neg_IHM_Fraction_zoom_mass.png')
