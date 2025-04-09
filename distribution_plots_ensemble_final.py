# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:53:06 2024

@author: Ji-Hoon Ha
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import random

# 데이터 파일 읽기
filename = '../goes_proton_data/0209_proton_1day.txt'  # 파일 이름을 적절히 변경
data = pd.read_csv(filename, delim_whitespace=True, comment='#',
                   names=["Year", "Month", "Day", "Time", "Julian_Day", "Seconds",
                          "P>1", "P>5", "P>10", "P>30", "P>50", "P>100"])

## 'UTC Date'와 'TIME'을 조합하여 시계열로 변환
#data['Datetime'] = pd.to_datetime(data[['Year', 'Month', 'Day', 'Time']].astype(str).agg(' '.join, axis=1), format='%Y %m %d %H%M')

t_axis = np.linspace(0,285*5,285)
t_axis = t_axis/60

spectrum_LB3 = np.load('spectrum_LB3.npy')

N_p = 100
delta_p = 0.1
gamma = 5. / 3.

c = 1.
m_i = 1
C_s = c * 1.e-4
p_th = m_i * C_s / np.sqrt(gamma)  # in the unit of m_i c
p_inj = 3 * p_th  # in the unit of p_th
grid_min = np.log10(p_th)

x_axis = np.zeros((N_p))
for j in range(N_p):
    x_axis[j] = grid_min + j * delta_p


accumulated_par = np.zeros((300))
accumulated_par_50 = np.zeros((300))
accumulated_par_100 = np.zeros((300))

accumulated_par_min = np.zeros((300))
accumulated_par_50_min = np.zeros((300))
accumulated_par_100_min  = np.zeros((300))

accumulated_par_max = np.zeros((300))
accumulated_par_50_max = np.zeros((300))
accumulated_par_100_max = np.zeros((300))

accumulated_par_std = np.zeros((300))
accumulated_par_50_std = np.zeros((300))
accumulated_par_100_std = np.zeros((300))

N_model = 100

delta_Ms=0.2

index = [0,10,20,50,70]
index = [0,10,15,20,35,50,60,70]

for ensb_i in range(N_model):
    # Applying random scaling factor
    random_number = random.uniform(-1, 1)
    perturbed_spectrum = spectrum_LB3.copy()
    
    for col in range(1, 8):  # Applying to columns 1 through 4
        perturbed_spectrum[:, index[col]] *= (10**x_axis)**(delta_Ms * random_number)
    
    # Save each perturbed spectrum with unique filename
    filename = f"spectrum_LB3_ensemble_{ensb_i}.npy"
    np.save(filename, perturbed_spectrum)
    print(f"Saved: {filename}")


# Define the color and label lists
colors = ['black','black','black','red', 'blue','green','magenta','magenta']  # 4 columns
labels = ['1.0h','2.0h','4.0hr','5.0hr', '6.0hr','7.0hr','8.0hr', '8.0hr']

# Initialize lists to store average, min, and max for each column
ensemble_avg_all = []
ensemble_std_all = []
ensemble_min_all = []
ensemble_max_all = []

# Iterate through each column (1, 2, 3, 4)
for col in range(1, 8):
    ensemble_data = []

    for ensb_i in range(N_model):
        # Load each ensemble model result
        filename = f"spectrum_LB3_ensemble_{ensb_i}.npy"
        spectrum_LB3 = np.load(filename)
        
        # Collect data for the current column scaled by the factor
        scaled_data = 4. * 3.14 * (10 ** (x_axis)) ** 2 * spectrum_LB3[:, index[col]] * (10 ** (x_axis)) ** ((-0.1))
        ensemble_data.append(scaled_data)

    # Convert to numpy array for easier calculations
    ensemble_data = np.array(ensemble_data)

    # Calculate the ensemble average, minimum, and maximum
    ensemble_avg = np.mean(ensemble_data, axis=0)
    ensemble_std = np.std(ensemble_data, axis=0)
    ensemble_min = np.min(ensemble_data, axis=0)
    ensemble_max = np.max(ensemble_data, axis=0)

    # Store the results for each column
    ensemble_avg_all.append(ensemble_avg)
    ensemble_std_all.append(ensemble_std)
    ensemble_min_all.append(ensemble_min)
    ensemble_max_all.append(ensemble_max)

    accumulated_par[col] = 0.7e+20*np.sum(ensemble_avg[52:])
    accumulated_par_50[col] = 0.7e+20*np.sum(ensemble_avg[52:])/10
    accumulated_par_100[col] = 0.7e+20*np.sum(ensemble_avg[52:])/100
    
    accumulated_par_min[col] = 0.7e+20*np.sum(ensemble_min[52:])
    accumulated_par_50_min[col] = 0.7e+20*np.sum(ensemble_min[52:])/10
    accumulated_par_100_min[col] = 0.7e+20*np.sum(ensemble_min[52:])/100
    
    accumulated_par_max[col] = 0.7e+20*np.sum(ensemble_max[52:])
    accumulated_par_50_max[col] = 0.7e+20*np.sum(ensemble_max[52:])/10
    accumulated_par_100_max[col] = 0.7e+20*np.sum(ensemble_max[52:])/100

    accumulated_par_std[col] = 0.7e+20*np.sum(ensemble_std[52:])
    accumulated_par_50_std[col] = 0.7e+20*np.sum(ensemble_std[52:])/10
    accumulated_par_100_std[col] = 0.7e+20*np.sum(ensemble_std[52:])/100
    

# Plot ensemble averages and error ranges (min and max) for all columns
plt.figure(figsize=(10, 6))

for i in range(2,7,1):
    plt.plot(x_axis, 0.7e+20*ensemble_avg_all[i], label=labels[i], color=colors[i])
    plt.fill_between(x_axis, 0.7e+20*(ensemble_min_all[i]), 0.7e+20*(ensemble_max_all[i]), color=colors[i], alpha=0.2)

plt.xlabel(r'$log~E [MeV]$', fontsize=15)
plt.ylabel(r'$EI(E)~{\rm[cm^{-2}s^{-1}sr^{-1}]}$', fontsize=15)
plt.yscale('log', base=10)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0, 3)
plt.ylim(1e-4, 1e+5)
plt.legend(fontsize=15)
#plt.show()
plt.savefig('fig_particle_dist_ensb_240209.png')

accumulated_par[1] = 0.0
accumulated_par_50[1] = 0.0
accumulated_par_100[1] = 0.0
accumulated_par_min[1] = 0.0
accumulated_par_50_min[1] = 0.0
accumulated_par_100_min[1] = 0.0
accumulated_par_max[1] = 0.0
accumulated_par_50_max[1] = 0.0
accumulated_par_100_max[1] = 0.0

accumulated_par[2] = 0.0
accumulated_par_50[2] = 0.0
accumulated_par_100[2] = 0.0
accumulated_par_min[2] = 0.0
accumulated_par_50_min[2] = 0.0
accumulated_par_100_min[2] = 0.0
accumulated_par_max[2] = 0.0
accumulated_par_50_max[2] = 0.0
accumulated_par_100_max[2] = 0.0

for idx in range(8):
    accumulated_par[idx] += 0.2
    accumulated_par_50[idx] += 0.2
    accumulated_par_100[idx] += 0.2

    accumulated_par_min[idx] += 0.2
    accumulated_par_50_min[idx] += 0.2
    accumulated_par_100_min[idx] += 0.2
    
    accumulated_par_max[idx] += 0.2
    accumulated_par_50_max[idx] += 0.2
    accumulated_par_100_max[idx] += 0.2

    accumulated_par_std[idx] += 0.2
    accumulated_par_50_std[idx] += 0.2
    accumulated_par_100_std[idx] += 0.2


    
for idx in range(8,300):
    accumulated_par[idx] = 0.2+accumulated_par[7]#*np.exp(-idx/tau_decay)
    accumulated_par_50[idx] =  0.2+accumulated_par_50[7]#*np.exp(-idx/tau_decay)
    accumulated_par_100[idx] =  0.2+accumulated_par_100[7]#*np.exp(-idx/tau_decay)

    accumulated_par_min[idx] = 0.2+accumulated_par_min[7]#*np.exp(-idx/tau_decay)
    accumulated_par_50_min[idx] =  0.2+accumulated_par_50_min[7]#*np.exp(-idx/tau_decay)
    accumulated_par_100_min[idx] =  0.2+accumulated_par_100_min[7]#*np.exp(-idx/tau_decay)
    
    accumulated_par_max[idx] = 0.2+accumulated_par_max[7]#*np.exp(-idx/tau_decay)
    accumulated_par_50_max[idx] =  0.2+accumulated_par_50_max[7]#*np.exp(-idx/tau_decay)
    accumulated_par_100_max[idx] =  0.2+accumulated_par_100_max[7]#*np.exp(-idx/tau_decay)

    accumulated_par_std[idx] = 0.2+accumulated_par_std[7]#*np.exp(-idx/tau_decay)
    accumulated_par_50_std[idx] =  0.2+accumulated_par_50_std[7]#*np.exp(-idx/tau_decay)
    accumulated_par_100_std[idx] =  0.2+accumulated_par_100_std[7]#*np.exp(-idx/tau_decay)

print(accumulated_par[4])

time = np.arange(0,300,1)

## 그림 사이즈와 해상도 설정
plt.figure(figsize=(10, 6), dpi=600)

# 라인 플롯 그리기
plt.plot(time,accumulated_par,linestyle='solid',marker='s',color='black',label='Model(>10MeV)')
plt.fill_between(time, accumulated_par_min, accumulated_par_max, color='black', alpha=0.2)
plt.plot(t_axis-12, data['P>10'][:285],color='red', label='Observation(>10MeV)')
#plt.plot(time,accumulated_par_50,linestyle='dashed',marker='^',color='blue',label='50MeV')
#plt.fill_between(time, accumulated_par_50_min, accumulated_par_50_max, color='blue', alpha=0.2)
#plt.plot(time,accumulated_par_100,linestyle='dotted',marker='o',color='green',label='100MeV')
#plt.fill_between(time, accumulated_par_100_min, accumulated_par_100_max, color='green', alpha=0.2)
# 플롯에 제목, 축 레이블 추가
#plt.title('Sample Line Plot')
plt.xlabel(r'$time [hour]$', fontsize=15)
plt.ylabel(r'$\Phi_{\rm obs}(E>E_0)~{\rm [cm^{-2}~s^{-1}~sr^{-1}]}$', fontsize=15)
plt.yscale('log',base=10)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0,10)
plt.ylim(1e-2,1e+4)
plt.legend(fontsize=15)

#plt.show()
plt.savefig('fig_particle_evol_ensb_240209.png')