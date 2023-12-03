# UNVERIFIED 
# filter T90 = n/a !!
# <= 6 implemented
# Verified

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


print("SWIFT")

file1 = open("Swift_GRB_table1.txt","r")
file1.readline()

i = 2

ra = []
dec = []
err = []
T90 = []
num = ['.','1','2','3','4','5','6','7','8','9','0']

while(i<1303): 
    a = file1.readline().split()
    if (len(a) == 5) & (a[1] != 'n/a') & (a[2] != 'n/a') & (a[3]!='n/a') & (a[4]!='n/a'):
        #if (a[4][0]) not in num:
            #rint(a[4][0])
            #if (float(a[4][1:])>=2):
            #    ra.append(float(a[1]))
            #    dec.append(float(a[2]))
        #else:
            #if (float(a[4])>=2):
        ra.append(float(a[1]))
        dec.append(float(a[2]))
        # T90.append(float(a[4]))
        if a[3][0] in num:
            err.append(float(a[3]))
        elif a[3][0] == 'n':    
            err.append(-1)
        else:
            err.append(float(a[3][1:]))

        if a[4][0] in num:
            T90.append(float(a[4]))
        elif a[4][0] == 'n':
            T90.append(-1)
        else:
            T90.append(float(a[4][1:]))
    i+=1
file1.close()

print(len(err))
print(len(T90))

################                <6 implementation
count = []
count1 = 0
for i in range(len(err)):
    if err[i]>6:
        count.append(i)

for i in count:
    ra.pop(i)
    dec.pop(i)
    err.pop(i)
    T90.pop(i)

print(len(ra))
print(len(T90))
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



plt.hist(L_err+S_err)
plt.title("Swift GRB Positional Uncertainity Histogram")
plt.xlabel("Error in degrees")
plt.ylabel("No. of GRBs")
plt.show()