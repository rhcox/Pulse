import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, arange
import sys


def getMax(data):
        max=data[0]
#        print 'max=', max
        for k in range(0,data.size):
#            print data[k]
            if data[k]>max:
                max=data[k]
#                print 'new max= ', max
#        print 'Final max=', max
        return max

def getMin(data):
        min=data[0]
#        print 'min=', min
        for k in range(0,data.size):
#            print data[k]
            if data[k]<min:
                min=data[k]
#                print 'new min= ', min
#        print 'Final min=', min
        return min
    
def getSlope(t0,t1,y0,y1):
    deltatime=t1-t0
    slope=0
    if deltatime==0:
        slope=float('nan')
    else:
        slope= (y1-y0)/deltatime
        #print 'Time= ',t1,'    slope= ', slope
    return slope

def FindNotch(xvector, yvector, start,end=-1):
    Peaks=[]
    minLen=min(len(xvector), len(yvector))
    if start<0:
        start=0
    if end> minLen or end==-1:
        end=minLen
    #print "xvector= ",xvector
    #print "yvector= ",yvector
    #print 'minLen= ',minLen
    
    slope=np.zeros(minLen)
    Notch=False
    slope[0]=getSlope(xvector[start],xvector[start+1],yvector[start],yvector[start+1])
    #print 'Slope 0= ',slope[0] 
    for i in range(start+1,minLen-1):
        slope[start-i]=getSlope(xvector[i],xvector[i+1],yvector[i],yvector[i+1])
        #print 'Slope ',i,'= ', slope[i]
        if slope[i-1]<0 and slope[i]>0 and Notch==False:
            #print 'Trough= ', i
            Peaks.append(i-1)
            Notch=True
        elif Notch==True:
            #print 'Looking for Peak'
            if slope[i-1]>0 and slope[i]<0:
                #print 'Peak=' , i
                Peaks.append(i-1)
                return Peaks
                break
    #print "Looking for Inflection"
    Peaks=[]
    i=start+1
    while i < minLen:
        #print 'slope= ',slope[i]
        if abs(slope[i])<0.02 :
            Peaks.append(i)
            #print "Append Slope"
            while abs(slope[i+1])<0.02 and i < minLen-2:
                #print 'slope= ',slope[i]
                i=i+1
            Peaks.append(i)
            #print 'Peaks= ',Peaks
            return Peaks    

    return Peaks
    
def inflection(xvector, yvector):
    Inflections=[]
    minLen=min(len(xvector), len(yvector))
    #print "xvector= ",xvector
    #print "yvector= ",yvector
    #print 'minLen= ',minLen
    
    slope=np.zeros(minLen)
    ZeroSlope=0
    for i in range(1,minLen):
        slope[i-1]=getSlope(xvector[i-1],xvector[i],yvector[i-1],yvector[i])
        #print 'Slope= ',slope[i-1]
        if slope[i-1]<0 and abs(slope[i-1])<0.06:
            ZeroSlope=ZeroSlope+1
        else:
            if ZeroSlope>=3:
                Inflections.append([i-ZeroSlope,i])
            ZeroSlope=0
    #print 'Inflction points: ', Inflections
    return Inflections


def PlotFTT(Freq,DataFFT,n):
        Y1 = DataFFT[range(n/2)]

        x1 = Freq
        y1 = abs(Y1)
        plt.scatter(x1, y1,label = "Frequency Distibution", color= "red", marker= "*", s=10)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Distibution [Arb. Units]')
        plt.title('PGP Frequency Domain')
       # plt.xlim([0,100])
       # plt.ylim([0,20])
        plt.show()
