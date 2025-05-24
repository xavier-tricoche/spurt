import numpy as np
import vtk
from vtk.util import numpy_support as nps
import argparse 
import sys
import nrrd
import time
import os
import vtk_camera
import vtk_colors
import regex as re

def main():
    parser = argparse.ArgumentParser(description='Visualize path of particles')
    parser.add_argument('-i', '--input', type=str, required=True, help='Nrrd namegroup')

    args = parser.parse_args()
    
    #Get all relevant files
    f = []
    #for (dirpath, dirnames, filenames) in os.walk("//home//ed3brute//spurt//spurt-main//build"):
    for (dirpath, dirnames, filenames) in os.walk("//mnt//c//Users//jacob//Desktop//lavd_results//reshaped - ORIGINAL"):
        for file in filenames:
            if(args.input in file and file.split('.')[1] == "nrrd"):
                f.append(file)
    #Get all relevant times
    t = []
    for i in f:
        time = i.split('_')[-1]
        time = time.split('h')[0]
        t.append(time)
    
    #Get all relevant sizes
    sizes = []
    for i in f:
        #print(i)
        data, header = nrrd.read(i)
        sizes.append(header['sizes'][2])
    
    for i in range(len(f)):
        #print("unu reshape -i " + f[i] + " -s " + str(sizes[i]) + " 3 -o //home//ed3brute//spurt//spurt-main//build//reshaped//reshaped_" + t[i] + "h.nrrd")
        #os.system("unu reshape -i " + f[i] + " -s " + str(sizes[i]) + " 3 -o //home//ed3brute//spurt//spurt-main//build//reshaped//reshaped_" + t[i] + "h.nrrd")
        print("unu reshape -i " + f[i] + " -s " + str(sizes[i]) + " 3 -o /mnt/c/Users/jacob/Desktop/lavd_results/reshaped - ORIGINAL/reshaped_" + t[i] + "h.nrrd")
        os.system("unu reshape -i " + f[i] + " -s " + str(sizes[i]) + " 3 -o \"/mnt/c/Users/jacob/Desktop/lavd_results/reshaped - ORIGINAL/reshaped_" + t[i] + "h.nrrd\"")
        
main()