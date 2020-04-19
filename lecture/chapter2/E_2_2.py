#coding=utf-8

from math import sqrt,log
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.diagnostic import acorr_ljungbox as lb_test

# Set the file path
filepath='data/m-ibm3dx2608.txt'
# Load the data
data=pd.read_csv(filepath,delimiter='\s+')
# Check the 1st row of the data
da=data.iloc[0,:]
# Get the IBM simple returns
sibm=data.iloc[:,2]
# Ljung-Box statistic Q(5) Box-Ljung test
sresult=lb_test(sibm,lags=5)
print(sresult)
# Log IBM returns
libm=np.log(sibm+1)
lresult=lb_test(libm,lags=5,boxpierce=True)
print(lresult)