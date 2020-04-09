import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
#import matplotlib.tri as tri
import os, subprocess

inputdir = '../gauss/vor_i_output_gnuplot/'
#data_path = '../gauss/vor_i_output_gnuplot/000002.dat'


npts=1216
count=0;
files = os.listdir(inputdir)
files.sort()

for file in files:
    data_path=inputdir+file
    print ('processing file :', data_path)

    x = []
    y = []
    z = []

    count=count+1
    if count > 2:
        break
        
# read data
    with open(data_path, 'r') as file:
        data = csv.reader(file, delimiter=' ')
        for row in data:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]))
        
    #print (x)

# ----------
# Tricontour
# ----------
# Directly supply the unordered, irregularly spaced coordinates
# to tricontour.
    fig, ax = plt.subplots()

    ax.tricontour(x, y, z, levels=20, linewidths=0.5, colors='k')
    #cntr = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

    #fig.colorbar(cntr, ax)
    #ax.plot(x, y, 'ko', ms=0.5)
    ax.set(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
    ax.set_title('tricontour (%d points)' % npts)

    output=str(count)+'.pdf'
    #plt.subplots_adjust(hspace=0.5)
    plt.show()
    print("outputing to "+output+"\n")
    plt.savefig(output)
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
#import matplotlib.tri as tri
import os, subprocess

inputdir = '../gauss/vor_i_output_gnuplot/'
#data_path = '../gauss/vor_i_output_gnuplot/000002.dat'


npts=1216
count=0
max=10
files = os.listdir(inputdir)
files.sort()

for file in files:
    data_path=inputdir+file
    print ('processing file :', data_path)

    x = []
    y = []
    z = []

    count=count+1
    if count > max:
        break
        
# read data
    with open(data_path, 'r') as file:
        data = csv.reader(file, delimiter=' ')
        for row in data:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]))
        
    #print (x)

# ----------
# Tricontour
# ----------
# Directly supply the unordered, irregularly spaced coordinates
# to tricontour.
    fig, ax = plt.subplots()

    ax.tricontour(x, y, z, levels=20, linewidths=0.5, colors='k')
    #cntr = ax.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

    #fig.colorbar(cntr, ax)
    #ax.plot(x, y, 'ko', ms=0.5)
    ax.set(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
    ax.set_title('tricontour (%d points)' % npts)

    output=str(count).zfill(6)+'.png'
    print("output to "+output+"\n")
    #fig=plt.figure()
    #plt.subplots_adjust(hspace=0.5)
    plt.show()
    fig.savefig(output)
