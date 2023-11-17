import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from random import sample, seed
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
    df.loc[:, 'IHM_Fraction'] = np.array(df['Intracluster_Stars_Mass'] / (0.17 * df['Mvir']))  # Modify this calculation as needed
    df.loc[:, 'IHM'] = np.array(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h)
    df.loc[:, 'hmass'] = np.array(df['Mvir'] * 1.0e10 / Hubble_h)
    df.loc[:, 'smass'] = np.array(df['Total_Stellar_Mass'] * 1.0e10 / Hubble_h)
    df.loc[:, 'log_IHM'] = np.array(np.log10(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h))
    df.loc[:, 'log_hmass'] = np.array(np.log10(df['Mvir'] * 1.0e10 / Hubble_h))
    df.loc[:, 'log_smass'] = np.array(np.log10(df['Total_Stellar_Mass'] * 1.0e10 / Hubble_h))
    df.loc[:, 'ColdGas'] = np.array(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h)
    df.loc[:, 'HotGas'] = np.array(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h)
    df.loc[:, 'EjectedMass'] = np.array(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h)
    df.loc[:, 'BlackHoleMass'] = np.array(df['Intracluster_Stars_Mass'] * 1.0e10 / Hubble_h)
    return df

# Function to save the plot with higher quality
def save_plot(save_filename):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, '/Users/michaelbradley/Documents/Honours/Images/', save_filename)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300)
    plt.show()

# -------------------------------------------------------------------------

def plot_IHM_vs_hmass(df, title, save_filename):
    plt.figure()  # New figure
    ax = plt.subplot(111)  # 1 plot on the figure

    ax.scatter('log_hmass', 'log_IHM')

    ax.set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$', fontsize=14)
    ax.set_ylabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$', fontsize=14)
    ax.set_xlim(11, 14.5)
    ax.set_ylim(6, 12.5)

    # Manually adding the title in the upper left corner
    ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
    # Create a custom legend with larger symbols
    leg = ax.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('large')

    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def BaryonFraction():
    
    print('Plotting the average baryon fraction vs halo mass')
    
    plt.figure()  # New figure
    ax = plt.subplot(111)  # 1 plot on the figure
    
    HaloMass = 'log_hmass'
    Baryons = 'smass' + 'ColdGas' + 'HotGas' + 'EjectedMass' + 'IHM' + 'BlackHoleMass'

    MinHalo = 11.0
    MaxHalo = 16.0
    Interval = 0.1
    Nbins = int((MaxHalo-MinHalo)/Interval)
    HaloRange = np.arange(MinHalo, MaxHalo, Interval)
    
    MeanCentralHaloMass = []
    MeanBaryonFraction = []
    MeanBaryonFractionU = []
    MeanBaryonFractionL = []

    MeanStars = []
    MeanCold = []
    MeanHot = []
    MeanEjected = []
    MeanICS = []
    MeanBH = []

    for i in range(Nbins-1):
        
        w1 = np.where((HaloMass >= HaloRange[i]) & (HaloMass < HaloRange[i+1]))[0]
        HalosFound = len(w1)
        
        if HalosFound > 2:  
            
            BaryonFraction = []
            CentralHaloMass = []
            
            Stars = []
            Cold = []
            Hot = []
            Ejected = []
            ICS = []
            BH = []
            
            for j in range(HalosFound):
                
                w2 = np.where(G.CentralGalaxyIndex == G.CentralGalaxyIndex[w1[j]])[0]
                CentralAndSatellitesFound = len(w2)
                
                if CentralAndSatellitesFound > 0:
                    BaryonFraction.append(sum(Baryons[w2]) / G.Mvir[w1[j]])
                    CentralHaloMass.append(np.log10(G.Mvir[w1[j]] * 1.0e10 / self.Hubble_h))

                    Stars.append(sum(G.StellarMass[w2]) / G.Mvir[w1[j]])
                    Cold.append(sum(G.ColdGas[w2]) / G.Mvir[w1[j]])
                    Hot.append(sum(G.HotGas[w2]) / G.Mvir[w1[j]])
                    Ejected.append(sum(G.EjectedMass[w2]) / G.Mvir[w1[j]])
                    ICS.append(sum(G.IntraClusterStars[w2]) / G.Mvir[w1[j]])
                    BH.append(sum(G.BlackHoleMass[w2]) / G.Mvir[w1[j]])                        
                            
            MeanCentralHaloMass.append(np.mean(CentralHaloMass))
            MeanBaryonFraction.append(np.mean(BaryonFraction))
            MeanBaryonFractionU.append(np.mean(BaryonFraction) + np.var(BaryonFraction))
            MeanBaryonFractionL.append(np.mean(BaryonFraction) - np.var(BaryonFraction))
            
            MeanStars.append(np.mean(Stars))
            MeanCold.append(np.mean(Cold))
            MeanHot.append(np.mean(Hot))
            MeanEjected.append(np.mean(Ejected))
            MeanICS.append(np.mean(ICS))
            MeanBH.append(np.mean(BH))
            
            print('  ', i, HaloRange[i], HalosFound, np.mean(BaryonFraction))
    
    plt.plot(MeanCentralHaloMass, MeanBaryonFraction, 'k-', label='TOTAL')#, color='purple', alpha=0.3)
    plt.fill_between(MeanCentralHaloMass, MeanBaryonFractionU, MeanBaryonFractionL, 
        facecolor='purple', alpha=0.25, label='TOTAL')
    
    plt.plot(MeanCentralHaloMass, MeanStars, 'k--', label='Stars')
    plt.plot(MeanCentralHaloMass, MeanCold, label='Cold', color='blue')
    plt.plot(MeanCentralHaloMass, MeanHot, label='Hot', color='red')
    plt.plot(MeanCentralHaloMass, MeanEjected, label='Ejected', color='green')
    plt.plot(MeanCentralHaloMass, MeanICS, label='ICS', color='yellow')
    # plt.plot(MeanCentralHaloMass, MeanBH, 'k:', label='BH')
    
    plt.xlabel(r'$\mathrm{Central}\ \log_{10} M_{\mathrm{vir}}\ (M_{\odot})$')  # Set the x-axis label
    plt.ylabel(r'$\mathrm{Baryon\ Fraction}$')  # Set the y-axis label
        
    # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        
    plt.axis([10.8, 15.0, 0.0, 0.23])

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# List of CSV files to process along with corresponding titles
csv_files = ['/Users/michaelbradley/Documents/Honours/millennium/tao.4489.0.csv']  # Add your CSV filenames here
titles = ['z = 0.00']

# Initialize lists to store DataFrames for old and new datasets
datasets = []

# Loop through the specified CSV files for 1 dataset
for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_dilute = dilute_dataframe(df_filtered, target_rows)
    df_calculated = perform_calculations(df_filtered)
    datasets.append(df_calculated)

print('Intrahalo mass vs. Halo mass')
plot_IHM_vs_hmass(datasets, titles, 'test.png')