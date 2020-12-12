#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 01:08:45 2020

@author: asif
"""

# import numpy as np
# import os
import subprocess

# os.system(". /home/asif/heasoft-6.28/x86_64-pc-linux-gnu-libc2.31/headas-init.sh")
# test = subprocess.Popen([". /home/asif/heasoft-6.28/x86_64-pc-linux-gnu-libc2.31/headas-init.sh"],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
p = subprocess.Popen(["bash"],stdin=subprocess.PIPE, stdout=subprocess.PIPE)
# P.stdin.write(b'cd\n')
# P.stdin.write(b'source .bashrc\n')
p.stdin.write(b'. /home/asif/heasoft-6.28/x86_64-pc-linux-gnu-libc2.31/headas-init.sh\n')
# output = P.communicate()[0]


no_of_specs=10
base_name="spec"
nhg=1.03e-2
nhz=1.093
z=1.238
for i in range(no_of_specs):
# i=0
    folder_name= base_name+str(i)
    data_name=folder_name+"pc.pi"
    
    
    
    with open(folder_name+'/'+folder_name+'_script.xcm','w') as f:
        f.write("log "+folder_name+".log\n")
        f.write("cd "+ folder_name+"\n")
        f.write("data "+ data_name+"\n")
        f.write("query yes \nsetplot energy \nignore bad \nstatistic cstat \nmodel phabs*zphabs*(powerlaw) & /*\n")
        f.write("newpar 1  "+ str(nhg)+" -1\n")
        f.write("newpar 2  "+ str(nhz)+" -1\n")
        f.write("newpar 3  "+ str(z)+" -1\n")
        f.write("renorm \nfit\nshow\n")
        f.write("err 4\n")
        f.write("err 5\n")
        f.write("save model "+folder_name+"_model\n")
        f.write("save data "+folder_name+"_data\n")
        f.write("exit")
    
    p.stdin.write(bytes("xspec<./"+folder_name+'/'+folder_name+'_script.xcm\n','utf-8'))

stdout, stderr = p.communicate()
print(stdout)


    
        
        
    