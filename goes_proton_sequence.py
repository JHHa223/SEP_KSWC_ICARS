# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:13:06 2024

@author: user
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 데이터 파일 읽기
filename = '0209_proton_1day.txt'  # 파일 이름을 적절히 변경
data = pd.read_csv(filename, delim_whitespace=True, comment='#',
                   names=["Year", "Month", "Day", "Time", "Julian_Day", "Seconds",
                          "P>1", "P>5", "P>10", "P>30", "P>50", "P>100"])

## 'UTC Date'와 'TIME'을 조합하여 시계열로 변환
#data['Datetime'] = pd.to_datetime(data[['Year', 'Month', 'Day', 'Time']].astype(str).agg(' '.join, axis=1), format='%Y %m %d %H%M')

t_axis = np.linspace(0,285*5,285)
t_axis = t_axis/60-3

# P>10 시계열 데이터 플로팅
plt.figure(figsize=(10, 6))
plt.plot(t_axis, data['P>10'][1:], label='10MeV')
plt.plot(t_axis, data['P>50'][1:], label='50MeV')
plt.plot(t_axis, data['P>100'][1:], label='100MeV')
plt.xlabel(r'$time [hour]$', fontsize=15)
plt.ylabel(r'${\rm particles~cm^{-2}~s^{-1}~sr^{-1}}$', fontsize=15)
plt.yscale('log',base=10)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(-3,12)
plt.ylim(1e-2,1e+4)
plt.legend(fontsize=15)
plt.show()
