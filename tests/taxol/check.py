#!/usr/bin/python

import sys
import re

infile = "taxol_rimp2_cc-pvdz.out"
#print infile

emp2_ref = -2921.33104029795

threshold = 0.0000001 

f = open(infile,'r')
lines = f.readlines()
f.close()

for line in lines:
    find = re.search("Total MP2 energy", line)
    if find:
        s1 = line.replace(' ','')
        s2 = s1.replace("TotalMP2energy=",'')
        s3 = s2.replace('\n','')
        emp2 = float(s3)

error = abs(emp2 - emp2_ref)
#print emp2, emp2_ref, error

if (error < threshold):
    print "The result is OK"
else:
    print "Error! The result is wrong!: error=", error
