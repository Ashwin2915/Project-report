from __future__ import print_function, division
# importing the csv library
import csv
import numpy as np
from matplotlib import pyplot as plt
from astroML.utils.decorators import pickle_results
from astroML.datasets import fetch_sdss_specgals
from astroML.correlation import bootstrap_two_point_angular
from astropy.io import fits
from scipy import stats
import random
import math


print("BATSE")

file1 = open("batse_grb.txt","r")
id1 = []
ra = []
dec = []
std_dev = []
stat = []

raf = []
decf = []
std_devf = []
statf = []
id1f = []  

for i in range(2702):
    a = file1.readline()
    ra.append(float(a.split()[5]))
    dec.append(float(a.split()[6]))
    std_dev.append(float(a.split()[9]))
    stat.append(a.split()[12])
    id1.append(float(a.split()[0]))
file1.close()

id2 = []
T90 = []

file2 = open("duration_table.txt","r")
for i in range(2041):
    a = file2.readline().split()
    id2.append(float(a[0]))
    T90.append(float(a[4]))

for i in range(len(std_dev)):
    if std_dev[i] <= 6:
        raf.append(ra[i])
        decf.append(dec[i])
        std_devf.append(std_dev[i])
        statf.append(stat[i])
        id1f.append(id1[i])


for i in range(len(statf)):
    if statf[i] == 'Y':
        raf.pop(i)
        decf.pop(i)
        std_devf.pop(i)
        id1f.pop(i)
id1 = []
T901 = []
L_ra = []
L_dec = []
L_err = []

S_ra = []
S_dec = []
S_err = []
j = 0

print("first cehck")
print(len(raf))
print(len(decf))

for i in range(len(id1f)):
    for j in range(len(id2)):
        if id1f[i] == id2[j]:
            id1.append(id1f[i])
            T901.append(T90[j])
            if T90[j] > 2:
                L_ra.append(raf[i])
                L_dec.append(decf[i])
                L_err.append(std_devf[i])
            else:
                S_ra.append(raf[i])
                S_dec.append(decf[i])
                S_err.append(std_devf[i])


raf = L_ra + S_ra
decf = L_dec + S_dec
std_devf = L_err + S_err 
print("check")
print(len(L_ra))
print(len(S_ra))
print(len(raf))
raf = L_ra
decf = L_dec
std_devf = L_err

bins = np.linspace(0, 180, 100)
result = []
for i in range(100):
    ra1 = []
    dec1 = []
    std1 = []
    for j in range(len(raf)):
        ra1.append(np.random.normal(raf[j],std_devf[j]))
        dec1.append(np.random.normal(decf[j],std_devf[j]))
    corr, dcorr, bootstraps = bootstrap_two_point_angular(ra1, dec1, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)
    result.append(corr)
    print(i)
corr = []
dcorr = []
corr_temp = []
for i in range(len(result[0])):
    for j in range(len(result)):
        corr_temp.append(result[j][i])
    corr.append(sum(corr_temp)/len(corr_temp))
    dcorr.append(np.std(corr_temp))
    corr_temp = []   
#corr, dcorr, bootstraps = bootstrap_two_point_angular(ra, dec, bins, method = 'landy-szalay', Nbootstraps = 10, random_state = None)








#BENCHMARK 
#BENCHMARK 
#BENCHMARK
rah = []
dech = []
for i in range(len(raf)):
    u = random.random()
    v = random.random()
    d = math.asin(2*u - 1) #changed
    dd = 90 - d*(180/math.pi)
    r = 360*v
    rah.append(r)
    dech.append(d)

bins = np.linspace(0, 180, 100)
corr_h, dcorr_h, bootstraps_h = bootstrap_two_point_angular(rah, dech, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)




absolute_sum_data = [0,0,0,0,0,0,0,0,0]
absolute_sum_homo = [0,0,0,0,0,0,0,0,0]
err_data = [0,0,0,0,0,0,0,0,0]
err_homo = [0,0,0,0,0,0,0,0,0]
sum_data = 0
sum_homo = 0
homo_err = 0
data_err = 0
x_axis = [1,2,3,4,5,6,7,8,9]


absolute_sum_data[0] = sum(corr[0:12])
absolute_sum_data[1] = sum(corr[12:23])
absolute_sum_data[2] = sum(corr[23:34])
absolute_sum_data[3] = sum(corr[34:45])
absolute_sum_data[4] = sum(corr[45:56])
absolute_sum_data[5] = sum(corr[56:67])
absolute_sum_data[6] = sum(corr[67:78])
absolute_sum_data[7] = sum(corr[78:89])
absolute_sum_data[8] = sum(corr[89:100])

absolute_sum_homo[0] = sum(corr_h[0:12])
absolute_sum_homo[1] = sum(corr_h[12:23])
absolute_sum_homo[2] = sum(corr_h[23:34])
absolute_sum_homo[3] = sum(corr_h[34:45])
absolute_sum_homo[4] = sum(corr_h[45:56])
absolute_sum_homo[5] = sum(corr_h[56:67])
absolute_sum_homo[6] = sum(corr_h[67:78])
absolute_sum_homo[7] = sum(corr_h[78:89])
absolute_sum_homo[8] = sum(corr_h[89:100])

err_data[0] = sum(dcorr[0:12])
err_data[1] = sum(dcorr[12:23])
err_data[2] = sum(dcorr[23:34])
err_data[3] = sum(dcorr[34:45])
err_data[4] = sum(dcorr[45:56])
err_data[5] = sum(dcorr[56:67])
err_data[6] = sum(dcorr[67:78])
err_data[7] = sum(dcorr[78:89])
err_data[8] = sum(dcorr[89:100])

err_homo[0] = sum(dcorr_h[0:12])
err_homo[1] = sum(dcorr_h[12:23])
err_homo[2] = sum(dcorr_h[23:34])
err_homo[3] = sum(dcorr_h[34:45])
err_homo[4] = sum(dcorr_h[45:56])
err_homo[5] = sum(dcorr_h[56:67])
err_homo[6] = sum(dcorr_h[67:78])
err_homo[7] = sum(dcorr_h[78:89])
err_homo[8] = sum(dcorr_h[89:100])

data_low=[]
data_high=[]
homo_low=[]
homo_high=[]

for i in range(len(absolute_sum_data)):
    data_low.append(absolute_sum_data[i] - err_data[i])
    data_high.append(absolute_sum_data[i] + err_data[i])
    homo_low.append(absolute_sum_homo[i] - err_homo[i])
    homo_high.append(absolute_sum_homo[i] + err_homo[i])

print("BATSE")
print("KS TEST")
KS_res = stats.ks_2samp(corr,corr_h)
if(KS_res.pvalue>0.05):
    print("KS test success")
else:
    print("KS test failed")
print(KS_res)

#AD test
print("AD test")
AD_res = stats.anderson_ksamp([corr,corr_h])
if(AD_res.pvalue>0.05):
    print("AD test success")
else:
    print("AD test failed")    
print(AD_res)

print("plotting")

#plotting the graph
plt.figure(figsize =(10,3))
plt.errorbar(x_axis, absolute_sum_data, yerr=err_data, xerr=None, fmt='ro', ecolor='red', elinewidth=None, capsize=10, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None,  data=None)
plt.errorbar(x_axis, absolute_sum_homo, yerr=err_homo, xerr=None, fmt='o', ecolor='blue', elinewidth=None, capsize=10, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None,  data=None)
plt.title("BATSE LGRB ABSOLUTE SUM TEST")
plt.xlabel("bins")
plt.ylabel("absolute sum")
plt.show()

