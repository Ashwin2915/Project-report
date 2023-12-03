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

# set ra = S_ra, dec = S_dec and err = S_err for SGRB analysis
ra = L_ra
dec = L_dec
err = L_err

raf = []
decf = []
errf = []


for i in range(len(err)):
    if err[i] <= 6:
        raf.append(ra[i])
        decf.append(dec[i])
        errf.append(err[i])
print("obtained filteres results")
#filtered values from the file <=6

bins = np.linspace(0, 180, 100)
result = []
for i in range(30):           
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


#########################      BENCHMARK
rah = []
dech = []
for i in range(len(ra)):  
    u = random.random()  
    v = random.random()
    #d = math.acos(1-2*u)
    #dd = 90 - d*(180/math.pi)
    d = math.asin(2*u - 1)
    r = 360*v
    rah.append(r)
    dech.append(d)
print("homo done")
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


print("plotting")

#plotting the graph
plt.figure(figsize =(10,3))
plt.errorbar(x_axis, absolute_sum_data, yerr=err_data, xerr=None, fmt='ro', ecolor='red', elinewidth=None, capsize=10, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None,  data=None)
plt.errorbar(x_axis, absolute_sum_homo, yerr=err_homo, xerr=None, fmt='o', ecolor='blue', elinewidth=None, capsize=10, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None,  data=None)
plt.title("FERMI LGRB ABSOLUTE SUM TEST")
plt.xlabel("bins")
plt.ylabel("absolute sum")
plt.show()