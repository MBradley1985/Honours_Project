#1 /usr/bin/env python
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

np.seterr(divide='ignore')

# Define the columns you want to extract and the target number of rows
columns_to_extract = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification']  # Change to your desired column names
columns_to_extract_2 = ['Mvir', 'Intracluster_Stars_Mass', 'Total_Stellar_Mass', 'Galaxy_Classification', 'Metals_IntraCluster_Stars_Mass']  # Change to your desired column names
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
    df.loc[:, 'halo_mass'] = df['Mvir'] * 1.0e10 / Hubble_h
    df.loc[:, 'IHM_Fraction'] = df['Intracluster_Stars_Mass'] / (0.17 * df['Mvir'])  # Modify this calculation as needed
    return df

def perform_calculations_2(df):

    df = df.copy ()
    df.loc[:, 'IHM_Fraction'] = df['Intracluster_Stars_Mass'] / (0.17 * df['Mvir'])  # Modify this calculation as needed
    df.loc[:, 'metallicity'] = np.log10((df['Metals_IntraCluster_Stars_Mass'] / df['Intracluster_Stars_Mass']) / 0.02) + 9.0
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

# Function to split up IHS into halo bins for mass function
def divide_and_store_custom_mass_ranges(df, bin_edges):
    adjusted_mvir = np.log10(df['Mvir'] * 1.0e10 / 0.73)
    mass_ranges = pd.cut(adjusted_mvir, bins=bin_edges)
    grouped = df.groupby(mass_ranges)['Intracluster_Stars_Mass'].apply(list)
    return grouped

# -------------------------------------------------------------------------

def Scatter_Plot(df, property_1, property_2, titles, ax):
    
    for df, title in zip(df, titles):
        # df_dilute = dilute_dataframe(df, target_rows)
        ax.scatter(np.log10(df[property_1] * 1.0e10 / Hubble_h), np.log10(df[property_2] * 1.0e10 / Hubble_h), alpha=0.6, s=0.8, color= 'gray')

# -------------------------------------------------------------------------

def calculate_statistics(df, property_to_bin, property_to_calculate):
    
    halo = np.log10(np.array(df[property_to_bin]) * 1.0e10 / Hubble_h)
    ihstars = np.log10(np.array(df[property_to_calculate]) * 1.0e10 / Hubble_h)

    mi = 6.5
    ma = 16.5
    bin_edges = np.arange(mi, ma + binwidth, binwidth)
                   
    mean_values = []
    std_values = []
    median_values = []
    mid_bin_values = []

    for i in range(len(bin_edges)):

       min_bin = mi + binwidth*i
       max_bin = mi + binwidth*(i+1)
       mid_bin_values.append( mi + binwidth*(i+0.5) )

       w = np.where((halo >= min_bin) & (halo < max_bin) & (ihstars > 0))[0] 
      
       mean_bin = np.mean(ihstars[w])
       mean_values.append(mean_bin)
       
       std_bin = np.std(ihstars[w])
       std_values.append(std_bin)
       
       median_bin = np.median(ihstars[w])
       median_values.append(median_bin)
       
       # print(i, min_bin, max_bin, mean_bin, std_bin) 
        
    
    return mean_values, std_values, median_values, mid_bin_values

# -------------------------------------------------------------------------


def IHS_MassFunction(data_per_file, datasets, property_1, subset_titles, main_titles, rows, cols, save_filename):
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
        ax.text(0.05, 0.95, main_title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        ax.set_ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$')

        handles_1 = custom_lines[:3]  # Assuming the first three lines correspond to your custom legend labels
        
        handles_2 = [Line2D([], [], color='black', label='Overall')]  # Line for the second histogram dataset
        
        ax.legend(handles=handles_1 + handles_2, loc='upper right', frameon=False, fontsize='x-small')

    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def Metallicity(df, property_1, property_2, property_3, titles, save_filename):
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
        Z = np.log10((df[property_2] / df[property_3]) / 0.02) + 9.0

        # overplot Tremonti et al. 2003 (h=0.7)
        w = np.arange(7.0, 13.0, 0.1)
        Zobs = -1.492 + 1.847*w - 0.08026*w*w

        if 'z = 0.00' in title:
            Tremonti_x = [8.609833118200260,8.711310621427260,8.815325062234930,8.924413378203950,
                          9.038575569334320,9.157811635626040,9.282121577079110,9.411505393693530,
                          9.543426147888620,9.675346902083720,9.807267656278810,9.939188410473910,
                          10.071109164669000,10.203029918864100,10.334950673059200,10.466871427254300,
                          10.598792181449400,10.730712935644500,10.862633689839600,10.994554444034700,
                          11.126475198229800,11.215268013553400]
            Tremonti_y = [8.453893713634990,8.500077060961780,8.545875777364410,8.591929890811790,
                          8.638012677468800,8.684156867976300,8.72975155462255,8.774586739271690,
                          8.817390979855460,8.857467009114060,8.894712519122800,8.929127509881690,
                          8.960848391956970,8.989602344216140,9.015696290433270,9.038993820042110,
                          9.059392625117970,9.076995013585530,9.091835088086360,9.103878745978900,
                          9.11292137141375,9.117603788971350]
            # Conversion from Kroupa IMF to Slapeter IMF to Chabrier IMF
            ax.plot(np.log10((10**w *1.5 /1.8)), Zobs, 'b', lw=2.0, label='Tremonti et al. 2004')
            # ax.scatter(Tremonti_x, Tremonti_y, c='r', label='Tremonti et al. 2004')

        ax.scatter(mass, Z, marker='o', s=1, c='gray', alpha=0.5, label='')
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

    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def IHMFraction_vs_redshift(df_list, property_1, property_2, titles, save_filename):
    plt.figure()
    
    colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
    redshifts = [0.00, 0.0199, 0.0414, 0.0645, 0.0893, 0.1159, 0.1749, 0.2075, 0.2798, 0.3197, 
                 0.4079, 0.5086, 0.5642, 0.6235, 0.6871, 0.7550, 0.827, 0.905, 1.077, 1.1734, 1.2758, 1.3857,
                   1.503, 1.6302, 1.766, 1.9126, 2.07]
    mean_values = []
    std_values = []
    max_values = []

    for df, color, title in zip(df_list, colors, titles):
        
        ihs_frac = np.array(df[property_2])
        halo = np.array(np.log10(df[property_1]))
        w = np.where((halo >= 14) & (ihs_frac > 0))[0]
        hmass = halo[w]
        ihs = ihs_frac[w]
        
        mean = np.mean(ihs)
        mean_values.append(mean)
        std = np.std(ihs)
        std_values.append(std)
        ma = np.max(ihs)
        max_values.append(ma)

    mean = np.array(mean_values) * 100
    std = np.array(std_values)
    maxx = np.array(max_values) * 100

    plt.plot(redshifts, mean, c ='r', label = 'SAGE (Millenium) Mean', linestyle = '--')
    plt.plot(redshifts, maxx, c ='r', label = 'SAGE (Millenium) Maximum', linestyle = ':')

    # Pentagons - Burke, Collins, Stott and Hilton - ICL @ z=1, 2012, 70
    redshifts_1 = [0.9468354430379745, 0.8303797468354429, 0.7949367088607594, 0.8075949367088605,  1.2227848101265821]
    ihs_fraction_1 = [1.415094339622641, 2.594339622641499, 3.7735849056603783, 1.5330188679245182, 2.3584905660377373]
    # Circles (gray) - Montes and Trujillo - ICL at the Frontier, 2018, 37
    redshifts_2 = [0.5341772151898734, 0.5443037974683542, 0.36708860759493667, 0.39746835443037976, 0.3417721518987341,
                   0.30379746835443033, 0.04810126582278479]
    ihs_fraction_2 = [1.5330188679245182, 0, 1.0613207547169807, 1.5330188679245182, 2.7122641509433976,
                      3.3018867924528266, 8.60849056603773]
    # Squares - Burke, Hilton and Collins, ICL and CLASH, 2015, 36
    redshifts_3 = [0.4025316455696200,0.3873417721518990,0.39746835443038000,0.33924050632911400,0.3443037974683540,
                   0.3417721518987340,0.2911392405063290,0.2253164556962030,0.2177215189873420,0.21265822784810100,
                   0.19493670886075900,0.17721518987341800]
    ihs_fraction_3 = [2.594339622641500,2.7122641509434000,3.3018867924528300,5.542452830188670,6.014150943396220,
                      7.193396226415100,12.971698113207500,12.500000000000000,16.27358490566040,18.042452830188700,
                      16.863207547169800,23.113207547169800]
    # Crosses Furnell et al., Growth of ICL in XCS-HSC from 0.1<z<0.5, 2021, 71
    redshifts_4 = [0.1443037974683540,0.12658227848101300,0.12151898734177200,0.08101265822784800,0.2253164556962030,
                   0.21518987341772100,0.2556962025316460,0.3063291139240510,0.260759493670886,0.2936708860759490,0.3215189873417720,
                   0.3417721518987340,0.3721518987341770,0.3367088607594940,0.37721518987341800,0.3291139240506330,0.4962025316455700,
                   0.42531645569620200,0.10886075949367100]
    ihs_fraction_4 = [38.561320754717000,30.660377358490600,31.014150943396200,28.891509433962300,26.533018867924500,
                      23.58490566037740,28.5377358490566,29.716981132075500,32.54716981132080,27.476415094339600,27.594339622641500,
                      26.650943396226400,19.81132075471700,18.867924528301900,15.448113207547200,15.330188679245300,
                      11.320754716981100,9.669811320754720,31.603773584905700]
    
    # Lines (model) - Rudick, Mihos and McBride, Quantity of ICL, 2011, 35
    redshifts_5 = [0.004011349760203840,0.023056412555524700,0.05248969142102060,0.08192297028651650,0.12001309587715800,
                   0.1581032214678000,0.19696284454512100,0.23774621133914200,0.2758363369297840,0.31392646252042500,0.3546136421286110,
                   0.38880818669293700,0.42681174381632700,0.4662869648829930,0.5043770904736340,0.542467216064276,0.5805573416549180,
                   0.6186474672455600,0.6567375928362010,0.6948277184268430,0.7329178440174850,0.7710079696081270,0.8090980951987680,
                   0.84718822078941,0.8852783463800520,0.9233684719706940,0.9614585975613350,0.999548723151977,1.0341761100525600,
                   1.0711119894131800,1.1449837481344300,1.183073873725070,1.2211639993157100,1.2592541249063500,1.2817619263917300,1.1083123425692700]
    ihs_fraction_5 = [17.035633925701400,15.255522799283600,13.375138908022200,11.340768199977500,11.708336986061200,12.328609312577500,
                      12.848617042445400,12.722826835652300,11.208443436987300,11.164335182657300,11.105524176883900,11.044875327180100,
                      10.871015291362500,10.55417099775830,10.091034327292800,9.55438389961052,8.833949078886400,8.635461934401190,9.150058234918420,
                      8.554596801462770,7.782702350686930,7.51070144898496,7.701837217748510,7.643026211975110,7.326917055943100,6.790266628260850,
                      6.6579418652707000,6.7755638768175,7.223997795839650,6.9544473527115800,6.15069694047515,5.900750165938210,
                      5.966912547433290,6.047777680371710,6.011020801763340,5.882352941176470]
    ihs_fraction_5_upper = [18.85060690943040,16.297065711530500,14.447949671672000,12.730913348946100,12.848617042445400,14.202658768425400,
                            14.648887274731100,14.773860661999500,13.173361147737800,12.298284887725600,12.632976687942900,12.607305217168800,
                            12.561096569775500,11.887220461955300,11.170216283234600,9.698470863755350,9.95870956430263,10.123590419774500,
                            9.023247003719530,8.291417550626820,8.370812408420900,8.436974789915970,7.682233549157370,6.9652293704367000,7.362623738019800,
                            7.385728061716490,6.192967350874780,6.091885934701760,6.253616200578600,6.280571244891410,8.352941176470590,
                            7.176470588235300,6.82352941176471,6.588235294117650,9.411764705882360,17.764705882352900]
    ihs_fraction_5_lower = [13.477568076410900,12.12532335338510,9.511133305781340,9.355896755125310,9.714643890343030,10.286764705882400,
                            10.192903390864500,9.31178850079526,9.795509023281450,10.201304963117900,9.694060038322330,9.213280066124820,
                            8.508854908083460,8.385515159864250,7.781068711637670,7.282808801613050,7.296776415484230,7.709188593470180,
                            6.988753772746070,6.318308306929330,7.069618905684490,7.481295946098260,7.046094503375120,6.706460945033760,
                            6.363886836403720,6.011020801763340,6.746158373930800,6.43740059362046,5.5905221104835400,5.202369472379120,
                            5.161936905909910,5.166429413295380,5.404532304725170,5.512352481976410,14.235294117647100,11.294117647058800]
    
    # Down Triangles - Feldmeier et al., Deep CCD, 2004, 72
    redshifts_6 = [0.16202531645569600,0.16202531645569600,0.16202531645569600,0.18481012658227800]
    ihs_fraction_6 = [15.212264150943400,12.146226415094300,10.259433962264100,7.311320754716980]
    # Black Circles - Montes and Trujillo - ICL at the Frontier, 2018,37
    redshifts_7 = [0.30126582278481000,0.38987341772151900,0.3417721518987340,0.5367088607594940,
                   0.5367088607594940,0.36962025316455700,0.043037974683544300]
    ihs_fraction_7 = [7.665094339622630,8.60849056603773,13.089622641509400,6.603773584905650,
                      5.778301886792450,4.834905660377360,10.849056603773600]
    # Star - Ko and Jee, Existence of ICL at z = 1.24, 2018, 74
    redshifts_8 = [1.2379746835443037]
    ihs_fraction_8 = [9.905660377358487]
    # Black Diamond - Kluge et al., ICL and host Cluster, 2021, 10
    redshifts_9 = [0.030379746835442978]
    ihs_fraction_9 = [17.924528301886788]
    # Plus - Zibetti et al., IGS in z=0.25 clusters, 2005, 73
    redshifts_10 = [0.24303797468354427]
    ihs_fraction_10 = [10.849056603773576]
    # Triangle - Presotto et al., ICL in CLASH-VLT cluster MACS J1206.2-0947, 2014, 75
    redshifts_11 = [0.4354430379746834]
    ihs_fraction_11 = [12.264150943396224]
     # Black Triangle - Presotto et al., ICL in CLASH-VLT cluster MACS J1206.2-0947, 2014, 75
    redshifts_12 = [0.43291139240506316]
    ihs_fraction_12 = [5.542452830188672]
     # Black Side Triangle - Spavone et al., Fornax Deep Survey, 2020, 76
    redshifts_13 = [0]
    ihs_fraction_13 = [34.08018867924528]

    plt.scatter(redshifts_1, ihs_fraction_1, marker = 'p', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_2, ihs_fraction_2, marker = 'o', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_3, ihs_fraction_3, marker = 's', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_4, ihs_fraction_4, marker = 'X', color = 'gray', edgecolors='black', s = 100)
    plt.plot(redshifts_5, ihs_fraction_5, linestyle = '--', color = 'plum')
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
        t.set_fontsize('large')
      
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def IHM_hmass(df, property_1, property_2, titles, save_filename):
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
        
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$', fontsize=14)
        ax.set_ylabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$', fontsize=14)
        ax.set_xlim(11, 14.5)
        ax.set_ylim(6, 12.5)
        
        # Manually adding the title in the upper left corner
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        
        # Calculate statistics for each bin
        mean_values, std_values, median_values, mid_bin_values = calculate_statistics(df, property_1, property_2)
        
        ax.xaxis.set_minor_locator(ticker.FixedLocator(mid_bin_values))
        
        mean_values = np.array(mean_values) 
        std_values = np.array(std_values)
        median_values = np.array(median_values)
        
        ax.plot(mid_bin_values, mean_values, color='blue', label='Mean')
        ax.plot(mid_bin_values, median_values, color='green', label='Median', linestyle = '--')
        ax.fill_between(mid_bin_values, (mean_values-std_values), (mean_values+std_values), color='orange', alpha=0.3, label=r'1$\sigma$ Std. Dev.')

        # Create a custom legend with larger symbols
        leg = ax.legend(loc='lower right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('large')
    
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------

def IHM_smass(df, property_1, property_2, titles, save_filename):
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
        
        ax.set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$', fontsize=14)
        ax.set_ylabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$', fontsize=14)
        ax.set_xlim(7, 12)
        ax.set_ylim(6, 12.5)
        
        # Manually adding the title in the upper left corner
        ax.text(0.05, 0.95, title, transform=ax.transAxes, fontsize=14, va='top', ha='left')
        
        # Calculate statistics for each bin
        mean_values, std_values, median_values, mid_bin_values = calculate_statistics(df, property_1, property_2)
        
        ax.xaxis.set_minor_locator(ticker.FixedLocator(mid_bin_values))
        
        mean_values = np.array(mean_values) 
        std_values = np.array(std_values)
        median_values = np.array(median_values)
        
        ax.plot(mid_bin_values, mean_values, color='blue', label='Mean')
        ax.plot(mid_bin_values, median_values, color='green', label='Median', linestyle = '--')
        ax.fill_between(mid_bin_values, (mean_values-std_values), (mean_values+std_values), color='orange', alpha=0.3, label=r'1$\sigma$ Std. Dev.')

        # Create a custom legend with larger symbols
        leg = ax.legend(loc='lower right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize('large')
    
    plt.tight_layout()
    save_plot(save_filename)

def Metallicity_all(save_filename):

    csv_files = '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4455.0.csv'
    titles = 'z = 0.00'

    print('Processing simulation data...', titles)
    df = pd.read_csv(csv_files)
    df = df[columns_to_extract_2]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_diluted = dilute_dataframe(df_filtered, target_rows)
    df_calculated = perform_calculations_2(df_diluted)

    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    
    axs[0, 0].scatter(np.log10(df_calculated['hmass']),df_calculated['metallicity'],s=1,c='gray')
    # axs[0, 0].set_title('Plot 1: sin(x)')
    axs[0, 0].set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
    axs[0, 0].set_xlabel(r'$\log_{10} M_{\mathrm{halo}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[0, 1].scatter(np.log10(df_calculated['smass']),df_calculated['metallicity'],s=1,c='gray')
    # axs[0, 1].set_title('Plot 2: cos(x)')
    axs[0, 1].set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
    axs[0, 1].set_xlabel(r'$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[1, 0].scatter(np.log10(df_calculated['IHM']),df_calculated['metallicity'],s=1,c='gray')
    # axs[1, 0].set_title('Plot 3: tan(x)')
    axs[1, 0].set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
    axs[1, 0].set_xlabel(r'$\log_{10} M_{\mathrm{IHS}}\ (M_{\odot})$')  # and the x-axis labels
    
    axs[1, 1].scatter(df_calculated['IHM_Fraction'],df_calculated['metallicity'],s=1,c='gray')
    # axs[1, 1].set_title('Plot 4: exp(x)')
    axs[1, 1].set_ylabel(r'$12\ +\ \log_{10}[\mathrm{O/H}]$')  # Set the y...
    axs[1, 1].set_xlabel('Intrahalo Stars Fraction')  # and the x-axis labels
    
    plt.tight_layout()
    save_plot(save_filename)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------    
# -------------------------------------------------------------------------

# List of CSV files to process along with corresponding titles
csv_files = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4455.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4456.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4457.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4459.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4460.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4461.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4462.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4463.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4464.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims_metallicity/tao.4465.0.csv']  # Add your CSV filenames here
titles = ['z = 0.00', 'z = 0.508', 'z = 1.0', 'z = 2.07', 'z = 3.06', 'z = 4.17', 'z = 5.288', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles
titles_sub = ['z = 0.00', 'z = 0.508', 'z = 1.0', 'z = 2.07', 'z = 3.06', 'z = 4.17', 'z = 5.288', 'z = 6.19', 'z = 8.54', 'z = 10.07']  # Add corresponding titles

# List of CSV files to process along with corresponding titles
csv_files2 = ['/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4404.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4406.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4427.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4407.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4428.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4408.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4409.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4410.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4411.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4412.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4413.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4414.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4415.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4416.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4417.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4418.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4419.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4429.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4420.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4430.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4431.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4432.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4421.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4433.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4434.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4435.0.csv',
             '/Users/michaelbradley/Documents/Honours/TAO/Small_sims/tao.4422.0.csv']  # Add your CSV filenames here
titles2 = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 0.1159', 'z = 0.1749', 'z = 0.2075',
            'z = 0.2798', 'z = 0.3197', 'z = 0.4079', 'z = 0.5086', 'z = 0.5642', 'z = 0.6235',
              'z = 0.6871', 'z = 0.7550', 'z = 0.827', 'z = 0.905', 'z = 1.077', 'z = 1.1734', 'z = 1.2758', 'z = 1.3857',
                'z = 1.503', 'z = 1.6302', 'z = 1.766', 'z = 1.9126', 'z = 2.07']  # Add corresponding titles


# Initialize lists to store DataFrames for old and new datasets
datasets = []
datasets2 = []
grouped_data_per_file = []

custom_bin_edges = [10.5, 12, 13.5, 17]  # Modify these edges as needed

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract_2]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_diluted = dilute_dataframe(df_filtered, target_rows)
    df_calculated = perform_calculations_2(df_diluted)
    datasets.append(df_calculated)
    mass_range_data = divide_and_store_custom_mass_ranges(df_calculated, custom_bin_edges)
    mass_range_arrays = [np.array(data) for data in mass_range_data]
    grouped_data_per_file.append(mass_range_arrays)

for idx, filename in enumerate(csv_files2):
    print('Processing simulation data...', titles2[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    datasets2.append(df_calculated)

IHS_MassFunction(grouped_data_per_file, datasets, columns_to_extract[1], titles_sub, titles, 5, 2, 'IHMfunction.png')
Metallicity(datasets, 'Total_Stellar_Mass', 'Metals_IntraCluster_Stars_Mass', 'Intracluster_Stars_Mass', titles, 'Metallicity_smass.png')
IHM_hmass(datasets, columns_to_extract[0], columns_to_extract[1], titles, 'IHM_vs_hmass.png')
IHM_smass(datasets, columns_to_extract[2], columns_to_extract[1], titles, 'IHM_vs_smass.png')
IHMFraction_vs_redshift(datasets2, 'halo_mass', 'IHM_Fraction', titles2, 'IHMFraction_vs_redshift.png')
Metallicity_all('Metallicity_all.png')