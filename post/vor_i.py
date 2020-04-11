import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv
import matplotlib.tri as tri
import os, subprocess

#inputdir = '../gauss/vor_i_output_gnuplot/'
inputdir = '../gauss/'
#inputdir = os.getcwd()
print(inputdir)

nx=160
ny=160
tot=nx*ny

x=np.zeros((ny, nx))
y=np.zeros((ny, nx))
z=np.zeros((ny, nx))

npts = 2508
fcount = 0
fmax   = 70

files = os.listdir(inputdir)
files.sort()
#data_path = '000074.grid'


zmin=0
zmax=8.0e-3
lev = np.linspace(zmin, zmax, 21)
#print(lev)

for file in files:
    #print (file)
    if file.endswith('.grid'):
    
        fcount = fcount+1
        if fcount > fmax:
            break
    
        data_path=inputdir+'/'+file
        print ('processing file :', data_path)

        
        #read data
        line_count=0
        with open(data_path, 'r') as file:
            data = csv.reader(file, delimiter=' ')
            for row in data:
                ix = int(line_count/ny)
                iy = line_count%ny
                #print(line_count, ix, iy)
                x[ix, iy]=float(row[0])
                y[ix, iy]=float(row[1])
                z[ix, iy]=float(row[2])
                line_count += 1

        fig, ax = plt.subplots()
        #ax.plot(x, y, 'ko', ms=1) #show grid points
        ax.set(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
        cp = ax.contour(x, y, z, levels=lev)
        fig.colorbar(cp) # Add a colorbar to a plot
        #ax.clabel(cp, inline=1, fontsize=10)
        title="Contour with grid " + str(nx) + '*' + str(ny)
        ax.set_title(title)
        output=str(fcount-1).zfill(6)+'.png'
        print("output to "+output+"\n")
        fig.savefig(output)
