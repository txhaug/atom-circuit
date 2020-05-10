# -*- coding: utf-8 -*-
"""
Created on Sat May  9 12:23:48 2020

@author: Tobias Haug
"""

import numpy as np
import matplotlib.pyplot as plt

figurenumber=0

#this function reorderes density in order to accomdate of the internal datastrucutre used to simulate rings and Y-junctions, whihc involves mapping ring into a next-nearest neighbor ordering. See 
def reorderFunction(array,systemLength,leadlength,model,shift=0):
    reorderedarray=np.array(array)
    if(model==0):#ring-lead
        reorderedarray[0+leadlength]=array[0+leadlength]
        for i in range(0,int(systemLength/2)+1):
            #print(i,max(i*2-1,0))
            reorderedarray[i+leadlength]=array[max(i*2-1,0)+leadlength]
        for i in range(0,int(np.ceil(systemLength/2))-1):
            #print(-i+systemLength-1,i*2+2)
            reorderedarray[-i+leadlength+systemLength-1]=array[i*2+2+leadlength]
    elif(model==1):#y-junction
        for i in range(0,systemLength):
            if(i>0):
                reorderedarray[i+leadlength]=array[i*2+leadlength]
            reorderedarray[i+leadlength+systemLength]=array[i*2+1+leadlength]



    return np.roll(reorderedarray,shift,axis=0)


def plot2D(griddata,gridx,gridy,xlabelstring,ylabelstring,cmap="RdYlBu_r",shading="flat"): #"YlOrBr" "cubehelix_r" "RdYlBu_r"
    global figurenumber
    plt.figure(figurenumber,figsize=(6,5))

    data=np.array(griddata)
    x=np.array(gridx)
    
    gridyIn=np.array(gridy)
    #if(inverty==True):
    #    gridyIn=np.flipud(gridyIn)
    y=np.array(gridyIn)

    deltax=(x[-1]-x[0])/len(x)
    deltay=(y[-1]-y[0])/len(y)
    if(shading=="flat"):
        addnum=1
    else:
        addnum=0
    plotx=np.linspace(x[0]-deltax/2,x[-1]+deltax/2,num=len(x)+addnum)
    ploty=np.linspace(y[0]-deltay/2,y[-1]+deltay/2,num=len(y)+addnum)

    
    plt.pcolormesh(plotx, ploty, data, cmap=cmap,linewidth=0,rasterized=True,shading=shading,antialiased=False)

    plt.colorbar()
        
    plt.xlabel(xlabelstring)
    plt.ylabel(ylabelstring)
            
    #show()
    figurenumber=figurenumber+1


#dataset to load

# uncomment for inputfileTest
#"""
#in this case, you should observe two density waves, one foreward and one backward, propagating in the chain. They are caused by an initial potential offset in the chain, which is quenched at t=0

dataset="DensityL20N20l20J1g1f0U5u5V0v0p0P0_1i2d2s1T10t0_01m7M100e3c9r1A10"

#this is the template file, already provided. Use this to check your results
#dataset="TemplateDensityL20N20l20J1g1f0U5u5V0v0p0P0_1i2d2s1T10t0_01m7M100e3c9r1A10"

N=20 #length of one part of chain
leadlength=20

model=2 #0: ring-lead, 1: y-junction 2:chain
#"""


# uncomment for inputfileYjunctionTest
"""
#in this case, you should observe two density waves, one foreward and one backward, propagating in the system. Site 0 to 40 represents the source lead of the Y-junction, and site 40:60 and 60:80 the two drain leads. Observe a density wave, that is transmitted and reflected at the site=40.
# Observe that the reflection is negative. Try changing the K paramter in the file, for values K<0.5 the reflection is positive. 

dataset="DensityL20N40l40J1g1f0U5u5V0v0p0P0_3i2d1s1T20t0_01m0M120e3c9r1A20"

#this is the template file, already provided. Use this to check your results
#dataset="TemplateDensityL20N40l40J1g1f0U5u5V0v0p0P0_3i2d1s1T20t0_01m7M120e3c9r1A20"
N=20 #length of two final parts of y-junction
leadlength=40

model=1 #0: ring-lead, 1: y-junction 2:chain
"""


if(model==0):#ring
    totallength=N+leadlength*2
elif(model==1):
    totallength=2*N+leadlength
elif(model==2):
    totallength=N+leadlength
    
    
datalist=np.loadtxt(dataset+".dat", delimiter=',')
#total length of system, for chain its leadlength+N, for ring its leadlength*2+N, for Y-junction its leadlength+2*N


timearray=datalist[0]

density=reorderFunction(datalist[8:totallength+8],N,leadlength,model)#density
densitsSq=reorderFunction(datalist[totallength+8:totallength*2+8],N,leadlength,model)#density**2


y=np.arange(totallength)
plot2D(density,timearray,y,"time","site") #plot density in time, for inputfileTest should generate forward and backward propagating density wave


        


        
    
    




