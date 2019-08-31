#coding=utf-8

from math import sqrt
import pandas as pd
import numpy as np
import scipy.stats as stats

# Set the file path
filepath='data/d-ibm3dx7008.txt'
# Load the data
ibm=pd.read_csv(filepath,delimiter='\s+').iloc[:,1]
# Percentage simple returns
sibm=ibm*100

# compute the summary statistics
summary=sibm.describe()
print(summary)

# Alternatively, one can use individual commands as follows:
# 算术平均值
mean=sibm.mean()
print('mean=%.8f'%mean)
# 中位数
median=sibm.median()
print('median=%.8f'%median)
# 方差
var=sibm.var()
print('var=%.8f'%var)
# 标准差
std=sibm.std()
print('std(sqrt)=%.8f'%std)
# 偏度
skewness=sibm.skew()
print('skewness=%.8f'%skewness)
# 峰度
kurtosis=sibm.kurt()
print('kurtosis=%.8f'%kurtosis)

# Simple tests
s1=sibm.skew()
# Compute test statistic
t1=s1/sqrt(6.0/len(sibm))
print('t1=%.8f'%t1)

# Compute p-value
pv = 2*(1-stats.norm.cdf(t1))
print('pv=%.8f'%pv)

# Turn to log returns in percentages
libm=np.log(ibm+1)*100
libm=libm.values
# Test mean being zero,One sample t-test
t,p_value=stats.ttest_1samp(libm,0)
print('One sample t-test:')
print('t=%.4f, p-value=%.4f\n'%(t,p_value))
# The result shows that the hypothesis of zero expected return
# cannot be rejected as the 5% or 10% level

# Normality test
statistic_x_squared,p_value=stats.jarque_bera(libm)
print('Normality test:')
print('STATISTIC:X-squared=%.4f,P-VALUE:Asymptotic_p_Value=%f\n'%(statistic_x_squared,p_value))
# The result shows the normality for log-return is rejuected