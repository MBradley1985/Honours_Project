import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d

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

def plot_IHMFraction_vs_hmass_single(df_list, property_1, property_2, titles, save_filename):
    plt.figure()
    
    colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
    
    for df, color, title in zip(df_list, colors, titles):
        if 'z = 0.00' in title:
            color = 'black'  # Set color to black for files with title 'z = 0.00'
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))
        w = np.where((halo >= 14) & (ihs_frac > 0))[0]
        hmass = halo[w]
        ihs = ihs_frac[w]
        
        plt.scatter(hmass, ihs, alpha=0.3, s=10.0, color=color, label=title)

    plt.ylim(0,0.4)    
    plt.xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')
    plt.ylabel(r'$F_{\mathrm{IHS}} \left(\frac{M_{\mathrm{IHS}}}{M_{\mathrm{halo}}} \times f_{\mathrm{cosmic\ baryons}}\right)$')
        
    # Create a custom legend with larger symbols
    leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')
      
    plt.tight_layout()
    
def plot_IHMFraction_vs_redshift(df_list, property_1, property_2, titles, save_filename):
    plt.figure()
    
    colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
    redshifts = [0.00, 0.0199, 0.0645, 0.1159, 0.1749, 0.201, 0.2798, 0.3197, 
                 0.4078, 0.508, 0.827, 1.077, 1.503, 2.07]
    mean_values = []
    std_values = []
    max_values = []

    for df, color, title in zip(df_list, colors, titles):
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))
        w = np.where((halo >= 14) & (ihs_frac > 0))[0]
        hmass = halo[w]
        ihs = ihs_frac[w]
        # mean_values = []
        
        mean = np.mean(ihs)
        mean_values.append(mean)
        std = np.std(ihs)
        std_values.append(std)
        ma = np.max(ihs)
        max_values.append(ma)
        print(mean, std, ma)

    mean = np.array(mean_values) * 100
    std = np.array(std_values)
    maxx = np.array(max_values) * 100
    print(mean)
    print(maxx)
    print(redshifts)

    plt.scatter(redshifts, mean, c ='r', label = 'Millenium Mean', marker = '+', s = 100)
    plt.scatter(redshifts, maxx, c ='r', label = 'Millenium Maximum', marker = 'x', s = 100)
    # plt.plot(redshifts, mean, c ='r', label = '')
    # plt.plot(redshifts, maxx, c ='r', label = '')
    # plt.fill_between(redshifts, (mean-std), (mean+std), color='red', alpha=0.8)

    # Pentagons
    redshifts_1 = [0.9468354430379745, 0.8303797468354429, 0.7949367088607594, 0.8075949367088605,  1.2227848101265821]
    ihs_fraction_1 = [1.415094339622641, 2.594339622641499, 3.7735849056603783, 1.5330188679245182, 2.3584905660377373]
    # Circles (gray)
    redshifts_2 = [0.5341772151898734, 0.5443037974683542, 0.36708860759493667, 0.39746835443037976, 0.3417721518987341,
                   0.30379746835443033, 0.04810126582278479]
    ihs_fraction_2 = [1.5330188679245182, 0, 1.0613207547169807, 1.5330188679245182, 2.7122641509433976,
                      3.3018867924528266, 8.60849056603773]
    # Squares
    redshifts_3 = [0.4025316455696200,0.3873417721518990,0.39746835443038000,0.33924050632911400,0.3443037974683540,
                   0.3417721518987340,0.2911392405063290,0.2253164556962030,0.2177215189873420,0.21265822784810100,
                   0.19493670886075900,0.17721518987341800]
    ihs_fraction_3 = [2.594339622641500,2.7122641509434000,3.3018867924528300,5.542452830188670,6.014150943396220,
                      7.193396226415100,12.971698113207500,12.500000000000000,16.27358490566040,18.042452830188700,
                      16.863207547169800,23.113207547169800]
    # Crosses
    redshifts_4 = [0.1443037974683540,0.12658227848101300,0.12151898734177200,0.08101265822784800,0.2253164556962030,
                   0.21518987341772100,0.2556962025316460,0.3063291139240510,0.260759493670886,0.2936708860759490,0.3215189873417720,
                   0.3417721518987340,0.3721518987341770,0.3367088607594940,0.37721518987341800,0.3291139240506330,0.4962025316455700,
                   0.42531645569620200,0.10886075949367100]
    ihs_fraction_4 = [38.561320754717000,30.660377358490600,31.014150943396200,28.891509433962300,26.533018867924500,
                      23.58490566037740,28.5377358490566,29.716981132075500,32.54716981132080,27.476415094339600,27.594339622641500,
                      26.650943396226400,19.81132075471700,18.867924528301900,15.448113207547200,15.330188679245300,
                      11.320754716981100,9.669811320754720,31.603773584905700]
    # Lines
    redshifts_5 = [0.0041887765689384000,0.0041887765689384000,0.03232097155110230,0.03899149201079060,0.0638505744692564,
                   0.06961788159962060,0.11158001278951100,0.1483129628693150,0.1884111503078090,0.22344588385260700,
                   0.2632949930335280,0.300011857824682,0.3374917052236000,0.37657783179675700,0.41301832208293100,
                   0.4534918329232500,0.4738079180624310,0.4963986092465620,0.5120909050484690,0.5438628840205930,
                   0.5842950913286540,0.6086851724743400,0.6391643198055320,0.6734717658353270,0.6852511464464150,0.7180538868379890,
                   0.7423543710906230,0.7802893309222420,0.8212791957759790,0.8561923960287580,0.8891269137585800,0.9264607357780210,
                   0.9647437227640590,1.0030267097501000,1.0413096967361300,1.0795926837221700,1.1178756707082100,1.1561586576942500,
                   1.1944416446802800,1.2327246316663200,1.2692674828802700,0.24556962025316500,0.5645569620253160]
    ihs_fraction_5 = [18.07154606914890,16.28979105490860,13.553541970453500,15.958198930670300,12.94292121426360,11.54528672743460,
                      11.743294595908300,12.753896846912800,13.095066051646600,13.840430563868500,10.748816480659700,10.721990853254400,
                      10.70827370896970,10.048490405357100,10.63321917977420,11.111419193927700,9.238162151749810,9.963925641397920,
                      9.083649432747120,9.444010244060510,7.950733914123380,7.671796769577460,8.711805708983250,8.94375014870684,
                      7.6593266384095300,7.1135750155305,6.652987053157180,7.454703122426600,7.696506103558380,8.107250356728790,
                      7.138901099342530,6.585194885535120,6.386239610992160,6.695725593614550,6.7988875878220100,6.098859769985660,
                      5.671474365411880,5.546206229588540,5.40620066602127,5.671474365411880,5.696527992576550,11.910377358490600,8.844339622641510]
    # Upside Down Triangles
    redshifts_6 = [0.16202531645569600,0.16202531645569600,0.16202531645569600,0.18481012658227800]
    ihs_fraction_6 = [15.212264150943400,12.146226415094300,10.259433962264100,7.311320754716980]
    # Black Circles
    redshifts_7 = [0.30126582278481000,0.38987341772151900,0.3417721518987340,0.5367088607594940,
                   0.5367088607594940,0.36962025316455700,0.043037974683544300]
    ihs_fraction_7 = [7.665094339622630,8.60849056603773,13.089622641509400,6.603773584905650,
                      5.778301886792450,4.834905660377360,10.849056603773600]
    # Star
    redshifts_8 = [1.2379746835443037]
    ihs_fraction_8 = [9.905660377358487]
    # Black Diamond
    redshifts_9 = [0.030379746835442978]
    ihs_fraction_9 = [17.924528301886788]
    # Plus
    redshifts_10 = [0.24303797468354427]
    ihs_fraction_10 = [10.849056603773576]
    # Triangle
    redshifts_11 = [0.4354430379746834]
    ihs_fraction_11 = [12.264150943396224]
     # Black Triangle
    redshifts_12 = [0.43291139240506316]
    ihs_fraction_12 = [5.542452830188672]
     # Black Side Triangle
    redshifts_13 = [0]
    ihs_fraction_13 = [34.08018867924528]

    plt.scatter(redshifts_1, ihs_fraction_1, marker = 'p', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_2, ihs_fraction_2, marker = 'o', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_3, ihs_fraction_3, marker = 's', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_4, ihs_fraction_4, marker = 'X', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_5, ihs_fraction_5, marker = '|', color = 'plum')
    plt.scatter(redshifts_6, ihs_fraction_6, marker = 'v', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_7, ihs_fraction_7, marker = 'o', color = 'k', edgecolors='gray', s = 100)
    plt.scatter(redshifts_8, ihs_fraction_8, marker = '*', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_9, ihs_fraction_9, marker = 'd', color = 'k', edgecolors='gray', s = 100)
    plt.scatter(redshifts_10, ihs_fraction_10, marker = 'P', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_11, ihs_fraction_11, marker = '^', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_12, ihs_fraction_12, marker = '^', color = 'k', edgecolors='gray', s = 100)
    plt.scatter(redshifts_13, ihs_fraction_13, marker = '<', color = 'k', edgecolors='gray', s = 100)

    plt.xlim(0, 2.1)
    plt.xlabel('Redshift')
    plt.ylabel('Intrahalo Stars %')

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
csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/tao.4358.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4391.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4392.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4393.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4394.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4395.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4396.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4397.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4364.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4365.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4366.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4367.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4368.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/tao.4369.0.csv']  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.0199', 'z = 0.0645', 'z = 0.1159', 'z = 0.1749', 'z = 0.201',
            'z = 0.2798', 'z = 0.3197', 'z = 0.4078', 'z = 0.508', 'z = 0.827', 'z = 1.077',
              'z = 1.503', 'z = 2.07']  # Add corresponding titles

datasets = []

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    datasets.append(df_calculated)

print('Intrahalo mass fraction vs. Halo mass')
plot_IHMFraction_vs_redshift(datasets, 'halo_mass', 'IHM_Fraction', titles, 'IHMFraction_vs_redshift.png')
