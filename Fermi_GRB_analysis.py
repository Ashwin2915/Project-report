#VERIFIED
#KD and AS tests there

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
import math
import random



print("FERMI")

#function for getting ra and dec list
def angle(a):
    ra = float(a[2])
    dec = float(a[3])

    if dec<0:
        dec = np.abs(dec)
        dec1 = dec//10000
        dec = dec%10000
        dec1 += (dec//100)/60
        dec = (dec%100)
        dec1 += dec/(60*60)
        dec1 = dec1 * -1
    else:
        dec1 = dec//10000
        dec = dec%10000
        dec1 += (dec//100)/60
        dec = (dec%100)
        dec1 += dec/(60*60)

    ra1 = (ra//10000)*360/24
    ra = ra%10000
    ra1 += (ra//100)*360/(24*60)
    ra = ra%100
    ra1 += ra*360/(24*60*60) 

    return [ra1,dec1]


#getting values from files 
file1 = open("Fermi_csv.txt","r")
file1.readline()
i=0
ra = []
dec = []
err = []
T90 = []
while(i<3613):
    a = file1.readline().split(',')
    lst = angle(a)
    ra.append(lst[0])
    dec.append(lst[1])
    if a[6][:len(a)-1] == '\n':
        err.append(-1)
    else:
        err.append(float(a[6][:len(a)-1]))
    T90.append(float(a[4]))
    i+=1
file1.close()
print(len(err))
print(len(ra))
c=0
for i in range(len(err)-1):
    if err[i] == -1:
        print(i)
        ra.pop(i)
        dec.pop(i)
        err.pop(i)
        T90.pop(i)
#getting 3 columns from file

L_ra = []
L_dec = []
L_err = []

S_ra = []
S_dec = []
S_err = []

for i in range(len(T90)):
    if T90[i] <= 2:
        S_ra.append(ra[i])
        S_dec.append(dec[i])
        S_err.append(err[i])

    else:
        L_ra.append(ra[i])
        L_dec.append(dec[i])
        L_err.append(err[i])


ra = L_ra 
dec = L_dec 
err = L_err 

raf = []
decf = []
errf = []

print(len(ra))
for i in range(len(err)):
    if err[i] <= 6:
        raf.append(ra[i])
        decf.append(dec[i])
        errf.append(err[i])
print("obtained filteres results")
#filtered values from the file <=6

bins = np.linspace(0, 180, 100)
result = []
for i in range(100):
    ra1 = []
    dec1 = []
    std1 = []
    for j in range(len(raf)):
        ra1.append(np.random.normal(raf[j],errf[j]**2))
        dec1.append(np.random.normal(decf[j],errf[j]**2))
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
print("done")

#########################      BENCHMARK
rah = []
dech = []
for i in range(len(ra)):  #COULD BE LEN(RAF)
    u = random.random()  
    v = random.random()
    #d = math.acos(1-2*u)
    #dd = 90 - d*(180/math.pi)
    d = math.asin(2*u - 1)
    r = 360*v
    rah.append(r)
    dech.append(d)

bins = np.linspace(0, 180, 100)
corr_h, dcorr_h, bootstraps_h = bootstrap_two_point_angular(rah, dech, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)
#corr1, dcorr1, boostraps1 = bootstrap_two_point_angular(ra, dec, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)
#corr2, dcorr2, boostraps2 = bootstrap_two_point_angular(ra, dec, bins, method = 'landy-szalay', Nbootstraps = 100, random_state = None)




print("FERMI")
########################        KS test 
print("KS TEST")
KS_res = stats.ks_2samp(corr_h,corr)
if(KS_res.pvalue > 0.05):
    print("KS test success")
else:
    print("KS test failed")
print(KS_res)

########################        AD test
print("AD test")
AD_res = stats.anderson_ksamp([corr_h,corr])
if(AD_res.pvalue > 0.05):
    print("AD test success")
else:
    print("AD test failed")    
print(AD_res)



#plotting the graph
bin_centers = 0.5 * (bins[1:] + bins[:-1])
plt.figure(figsize =(10,3))
plt.plot(bin_centers, corr, 'o', markersize = '4', color = 'black')
plt.fill_between(bin_centers, corr_h+dcorr_h, corr_h-dcorr_h, color= 'skyblue', alpha = 0.5)
plt.title("FERMI LGRB <6 2pacf")
plt.xlabel("deg")
plt.ylabel("w(theta)")
plt.ylim(-0.5, 0.5)
plt.errorbar(bin_centers, corr, yerr=dcorr, xerr=None, fmt='none', ecolor='gray', elinewidth=1)
plt.show