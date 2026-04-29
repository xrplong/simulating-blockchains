


# Importing
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from functions import predict_speed



### THIS FILE PLOTS RESULTS ###





# averaging samples in data with same parameters
df = pd.read_csv("vN_speed_simulations.csv")
group_cols = df.columns[:-1]
averaged_df = df.groupby(list(group_cols), as_index=False)[df.columns[-1]].mean()

#print(averaged_df)

# Renaming column label N to 1/N and inverting values
averaged_df['N'] = 1 / averaged_df['N']
averaged_df = averaged_df.rename(columns={'N': '1/N'})

#print(averaged_df)


# PLOTTING

if 1:
    # compute Derrida speed prediction v0
    lam=1
    omega=1
    T=1000
    g0, v0 = predict_speed(lam, omega)
    print('')
    print(f'lambda={lam}, omega={omega}, T={T}')
    print('v0 = ', v0)
    print('g0 = ', g0)

    # Getting dataframe for the triple (lambda, omega, T)
    df = averaged_df[averaged_df["omega"] == omega]
    print(df)

    # compute Derrida leading order 1/log(N)^2 coefficient
    dd_v_g0 = ( lam*( -2 + math.exp(g0)*( 2 - 2*g0 + g0**2 ) ) + 2*omega ) / ( g0**3 )
    Derrida_logN2_coef = ( math.pi**2 * g0**2 ) * dd_v_g0 / 2

    # getting axis values
    inverted_N = df[['1/N']]
    averaged_v = df[['v_N']]
    leading_correction_v0 = pd.DataFrame({'v0-v_N': v0-averaged_v['v_N']})

    Derrida_logN2_coef = round(Derrida_logN2_coef,3)
    Derrida_logN2_prediction = pd.DataFrame({'(1/2)v''(g0)pi^2g0^2/log(N)^2': Derrida_logN2_coef / (np.log(inverted_N['1/N']) ** 2)})

    # plotting v0 - v_N against 1/N as loglog plot
    plt.loglog(inverted_N, leading_correction_v0, ls=':', lw='2', marker='*', ms='10', color='k', label='simulated')
    plt.loglog(inverted_N, Derrida_logN2_prediction, ls='-', color='k', label='prediction')

    plt.xlabel('1/N', fontsize=12)
    plt.ylabel('$v_0-v_{1/N}$', fontsize=12)
    plt.legend()
    plt.grid(True)
    plt.show()
