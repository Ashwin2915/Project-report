# Analysis of the Batse GRB
# Unverified 

# table length verified
#ra verified 
#dec verified
#std_dev verified

#std <= 6 implemented  
#NO REDSHIFT

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


print(len(id2))
print(len(id1f))
print(len(id1))
print(len(L_ra))
print(len(S_ra))


raf = L_ra + S_ra
decf = L_dec + S_dec
std_devf = L_err + S_err 
print("check")
print(len(raf))
print(len(decf))
print(len(std_devf))


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
for i in range(4000):
    u = random.random()
    v = random.random()
    d = math.asin(2*u - 1) #changed
    dd = 90 - d*(180/math.pi)
    r = 360*v
    rah.append(r)
    dech.append(d)

bins = np.linspace(0, 180, 100)
corr_h, dcorr_h, bootstraps_h = bootstrap_two_point_angular(rah, dech, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)







#KS test 
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







#PLOTTING GRAPH  
bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.figure(figsize =(10,3))
plt.plot(bin_centers, corr, 'o', markersize = '4', color = 'black')
plt.errorbar(bin_centers, corr, yerr=dcorr, xerr=None, fmt='none', ecolor='gray', elinewidth=1)
plt.fill_between(bin_centers, corr_h+dcorr_h, corr_h-dcorr_h, color= 'skyblue', alpha = 0.5)
plt.title("Batse_LGRB_Analysis <6 vs benchmark")
plt.xlabel("deg")
plt.ylabel("w(theta)")
plt.show()