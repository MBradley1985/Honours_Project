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
    # print(mean)
    # print(maxx)
    print(redshifts)

    plt.plot(redshifts, mean, c ='r', label = 'SAGE (Millenium) Mean', linestyle = '--')
    plt.plot(redshifts, maxx, c ='r', label = 'SAGE (Millenium) Maximum', linestyle = ':')

    # Pentagons - Burke, Collins, Stott and Hilton - ICL @ z=1, 2012
    redshifts_1 = [0.9468354430379745, 0.8303797468354429, 0.7949367088607594, 0.8075949367088605,  1.2227848101265821]
    ihs_fraction_1 = [1.415094339622641, 2.594339622641499, 3.7735849056603783, 1.5330188679245182, 2.3584905660377373]
    # Circles (gray) - Montes and Trujillo - ICL at the Frontier, 2018
    redshifts_2 = [0.5341772151898734, 0.5443037974683542, 0.36708860759493667, 0.39746835443037976, 0.3417721518987341,
                   0.30379746835443033, 0.04810126582278479]
    ihs_fraction_2 = [1.5330188679245182, 0, 1.0613207547169807, 1.5330188679245182, 2.7122641509433976,
                      3.3018867924528266, 8.60849056603773]
    # Squares - Burke, Hilton and Collins, ICL and CLASH, 2015
    redshifts_3 = [0.4025316455696200,0.3873417721518990,0.39746835443038000,0.33924050632911400,0.3443037974683540,
                   0.3417721518987340,0.2911392405063290,0.2253164556962030,0.2177215189873420,0.21265822784810100,
                   0.19493670886075900,0.17721518987341800]
    ihs_fraction_3 = [2.594339622641500,2.7122641509434000,3.3018867924528300,5.542452830188670,6.014150943396220,
                      7.193396226415100,12.971698113207500,12.500000000000000,16.27358490566040,18.042452830188700,
                      16.863207547169800,23.113207547169800]
    # Crosses Furnell et al., Growth of ICL in XCS-HSC from 0.1<z<0.5, 2021
    redshifts_4 = [0.1443037974683540,0.12658227848101300,0.12151898734177200,0.08101265822784800,0.2253164556962030,
                   0.21518987341772100,0.2556962025316460,0.3063291139240510,0.260759493670886,0.2936708860759490,0.3215189873417720,
                   0.3417721518987340,0.3721518987341770,0.3367088607594940,0.37721518987341800,0.3291139240506330,0.4962025316455700,
                   0.42531645569620200,0.10886075949367100]
    ihs_fraction_4 = [38.561320754717000,30.660377358490600,31.014150943396200,28.891509433962300,26.533018867924500,
                      23.58490566037740,28.5377358490566,29.716981132075500,32.54716981132080,27.476415094339600,27.594339622641500,
                      26.650943396226400,19.81132075471700,18.867924528301900,15.448113207547200,15.330188679245300,
                      11.320754716981100,9.669811320754720,31.603773584905700]
    
    # Lines (model) - Rudick, Mihos and McBride, Quantity of ICL, 2011
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
    
    # Down Triangles - Feldmeier et al., Deep CCD, 2004
    redshifts_6 = [0.16202531645569600,0.16202531645569600,0.16202531645569600,0.18481012658227800]
    ihs_fraction_6 = [15.212264150943400,12.146226415094300,10.259433962264100,7.311320754716980]
    # Black Circles - Montes and Trujillo - ICL at the Frontier, 2018
    redshifts_7 = [0.30126582278481000,0.38987341772151900,0.3417721518987340,0.5367088607594940,
                   0.5367088607594940,0.36962025316455700,0.043037974683544300]
    ihs_fraction_7 = [7.665094339622630,8.60849056603773,13.089622641509400,6.603773584905650,
                      5.778301886792450,4.834905660377360,10.849056603773600]
    # Star - Ko and Jee, Existence of ICL at z = 1.24, 2018
    redshifts_8 = [1.2379746835443037]
    ihs_fraction_8 = [9.905660377358487]
    # Black Diamond - Kluge et al., ICL and host Cluster, 2021
    redshifts_9 = [0.030379746835442978]
    ihs_fraction_9 = [17.924528301886788]
    # Plus - Zibetti et al., IGS in z=0.25 clusters, 2005
    redshifts_10 = [0.24303797468354427]
    ihs_fraction_10 = [10.849056603773576]
    # Triangle - Presotto et al., ICL in CLASH-VLT cluster MACS J1206.2-0947, 2014
    redshifts_11 = [0.4354430379746834]
    ihs_fraction_11 = [12.264150943396224]
     # Black Triangle - Presotto et al., ICL in CLASH-VLT cluster MACS J1206.2-0947, 2014
    redshifts_12 = [0.43291139240506316]
    ihs_fraction_12 = [5.542452830188672]
     # Black Side Triangle - Spavone et al., Fornax Deep Survey, 2020
    redshifts_13 = [0]
    ihs_fraction_13 = [34.08018867924528]

    plt.scatter(redshifts_1, ihs_fraction_1, marker = 'p', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_2, ihs_fraction_2, marker = 'o', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_3, ihs_fraction_3, marker = 's', color = 'gray', edgecolors='black', s = 100)
    plt.scatter(redshifts_4, ihs_fraction_4, marker = 'X', color = 'gray', edgecolors='black', s = 100)
    plt.plot(redshifts_5, ihs_fraction_5, linestyle = '--', color = 'plum')
    # plt.fill_between(redshifts_5, ihs_fraction_5_upper, ihs_fraction_5_lower, color='purple', alpha=0.3)
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

# def plot_IHMFraction_vs_redshift2(df_list, property_1, property_2, titles, save_filename):
#     plt.figure()
    
#     colors = plt.cm.RdYlBu_r(np.linspace(0, 1, len(df_list)))
#     redshifts = [0.00, 0.0199, 0.0414, 0.0645, 0.0893, 0.1159, 0.1749, 0.2075, 0.2798, 0.3197, 
#                  0.4079, 0.5086, 0.5642, 0.6235, 0.6871, 0.7550, 0.827, 0.905, 1.077, 1.1734, 1.2758, 1.3857,
#                    1.503, 1.6302, 1.766, 1.9126, 2.07]
#     mean_values = []
#     std_values = []
#     max_values = []

#     for df, color, title in zip(df_list, colors, titles):
        
#         ihs_frac = np.array(df[property_2])
#         halo = np.array(np.log10(df[property_1]))
#         w = np.where((halo >= 14) & (ihs_frac > 0))[0]
#         hmass = halo[w]
#         ihs = ihs_frac[w]
#         # mean_values = []
        
#         mean = np.mean(ihs)
#         mean_values.append(mean)
#         std = np.std(ihs)
#         std_values.append(std)
#         ma = np.max(ihs)
#         max_values.append(ma)
#         print(mean, std, ma)

#     mean = np.array(mean_values)
#     std = np.array(std_values)
#     maxx = np.array(max_values)
#     # print(mean)
#     # print(maxx)
#     print(redshifts)

#     plt.plot(redshifts, mean, c ='k', label = 'SAGE (Millenium) Mean', linestyle = '--')
#     plt.plot(redshifts, maxx, c ='k', label = 'SAGE (Millenium) Maximum', linestyle = ':')

#     redshifts_1  = [0.01399,0.04977,0.07049,0.10439,0.14771,0.16466,0.17973,0.20610,
#                     0.23906,0.28143,0.33280,0.37969,0.41328,0.43911,0.46978,0.51666,
#                     0.55014,0.56301,0.60265,0.65511,0.69956,0.73063,0.74853,0.78620,
#                     0.83328,0.88639,0.92557,0.96754,1.01033,1.04440,1.09508,1.14405,
#                     1.18549,1.22795,1.27542,1.31904,1.34558,1.35876,1.37006,1.38890,
#                     1.41937,1.46047,1.48307,1.51245,1.55841,1.59984,1.64881,1.68460,
#                     1.72038,1.77061,1.79949,1.80702,1.82020,1.83715,1.85222,1.86729,
#                     1.88801,1.91814,1.96146,1.99348]
#     ihs_1 = [0.07899,0.07444,0.06466,0.05709,0.05802,0.06983,0.07775,0.08507,0.09206,
#              0.09909,0.10192,0.10311,0.10730,0.11685,0.11836,0.11351,0.11809,0.12474,
#              0.13488,0.13154,0.12283,0.12956,0.14013,0.14427,0.14177,0.14862,0.15378,
#              0.15559,0.16523,0.17331,0.17428,0.17372,0.17548,0.17994,0.17402,0.16889,
#              0.18168,0.18848,0.19727,0.20716,0.21191,0.20989,0.20074,0.19464,0.19916,
#              0.20546,0.20534,0.20097,0.19207,0.18985,0.19232,0.20099,0.20940,0.21809,
#              0.22610,0.23543,0.24327,0.25209,0.25698,0.25827]
#     ihs_1.reverse()

#     redshifts_2 = [0.03094,0.08744,0.13829,0.17596,0.20107,0.24430,0.30592,0.36631,
#                    0.41328,0.46576,0.52440,0.58484,0.64431,0.70144,0.75418,0.80263,
#                    0.86449,0.94095,0.99149,1.04649,1.10814,1.16916,1.22953,1.28945,
#                    1.33239,1.35782,1.37948,1.40396,1.44163,1.49625,1.54145,1.60172,
#                    1.65069,1.70657,1.75843,1.79158,1.81832,1.83150,1.85787,1.88236,
#                    1.90308,1.93823,1.98595]
#     ihs_2 = [0.05451,0.04560,0.04302,0.05548,0.06509,0.07169,0.07676,0.07738,0.08409,
#              0.09197,0.08813,0.09781,0.10087,0.09786,0.10481,0.11316,0.11081,0.12250,
#              0.12615,0.13656,0.13816,0.14027,0.14122,0.13971,0.13560,0.14763,0.15616,
#              0.16608,0.17167,0.16581,0.15990,0.16337,0.16215,0.15066,0.14851,0.14518,
#              0.15822,0.16717,0.18076,0.19065,0.20066,0.21076,0.21516]
#     ihs_2.reverse()

#     plt.plot(redshifts_1, ihs_1, c ='r', label = '', linestyle = '--')
#     plt.plot(redshifts_2, ihs_2, c ='b', label = '', linestyle = '--')

#     # plt.xlim(0, 2.1)
#     plt.xlabel('Redshift')
#     plt.ylabel('Intrahalo Stars %')

#     # Create a custom legend with larger symbols
#     leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
#     leg.draw_frame(False)  # Don't want a box frame
#     for t in leg.get_texts():  # Reduce the size of the text
#         t.set_fontsize('medium')
      
#     plt.tight_layout()
#     save_plot(save_filename)


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
titles = ['z = 0.00', 'z = 0.0199', 'z = 0.0414', 'z = 0.0645', 'z = 0.0893', 'z = 0.1159', 'z = 0.1749', 'z = 0.2075',
            'z = 0.2798', 'z = 0.3197', 'z = 0.4079', 'z = 0.5086', 'z = 0.5642', 'z = 0.6235',
              'z = 0.6871', 'z = 0.7550', 'z = 0.827', 'z = 0.905', 'z = 1.077', 'z = 1.1734', 'z = 1.2758', 'z = 1.3857',
                'z = 1.503', 'z = 1.6302', 'z = 1.766', 'z = 1.9126', 'z = 2.07']  # Add corresponding titles

datasets = []

for idx, filename in enumerate(csv_files):
    print('Processing simulation data...', titles[idx])
    df = pd.read_csv(filename)
    df = df[columns_to_extract]
    df_filtered = df[df['Galaxy_Classification'] == 0]
    df_calculated = perform_calculations(df_filtered)
    datasets.append(df_calculated)

print('Intrahalo mass fraction vs. Redshift')
plot_IHMFraction_vs_redshift(datasets, 'halo_mass', 'IHM_Fraction', titles, 'IHMFraction_vs_redshift.png')

# print('Intrahalo mass fraction vs. Redshift')
# plot_IHMFraction_vs_redshift2(datasets, 'halo_mass', 'IHM_Fraction', titles, 'IHMFraction_vs_redshift2.png')