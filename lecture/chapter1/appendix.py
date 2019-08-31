#coding=utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Set the file path
filepath='data/m-gm3dx7508.txt'
# Load the data
df=pd.read_csv(filepath,delimiter='\s+',parse_dates=['date'])
# Without Time Series in Column 2 contains GM stock returns
gm=df['gm']
# Set the Time Series
gm1=df.set_index('date')
# Column 2 contains GM stock returns
gm1=gm1['gm']

# Put four plots on a page
plt.figure()

ax1 = plt.subplot2grid((4,1),(0,0))
gm.plot(ax=ax1)
ax2 = plt.subplot2grid((4,1),(1,0))
gm1.plot(ax=ax2)

ax3 = plt.subplot2grid((4,1),(2,0))
sm.graphics.tsa.plot_acf(gm,lags=24,ax=ax3)
ax4 = plt.subplot2grid((4,1),(3,0))
sm.graphics.tsa.plot_acf(gm1,lags=24,ax=ax4)

plt.show()