import time
import math
from scipy.signal import find_peaks_cwt
import numpy as np
#from pypeaks import Data, Intervals
import csv
#import plotly
#import plotly.plotly as py
#import plotly.graph_objs as go
import matplotlib.pyplot as plt
from numpy import sin, linspace, pi
from scipy import fft, arange
import sys
import PPGfunc

class PPGWave(object):
    def ISO2mSec(self, rData):
        dateString=rData[4]
        #print 'Time '+dateString
        #print 'hour '+dateString[11:13]
        hour=int(dateString[11:13])
        #print 'Minute '+dateString[14:16]
        minute=int(dateString[14:16])
        #print 'Second '+dateString[17:19]
        second=int(dateString[17:19])
        #print 'milliSecond '+dateString[20:23]
        msecond=int(dateString[20:23])
        return (hour*12*60*60+minute*60+second)*1000+msecond
        
    def make(self, openF):
        with open(openF) as f:
            ppgData = csv.reader(f)
            next(ppgData) # skip header
            self.rawData = [raw for raw in ppgData]
        self.radius=3
        self.length=len(self.rawData)
        self.Time=np.zeros(self.length)
        self.Data=np.zeros(self.length)
        self.DataSN=np.zeros(self.length)
        self.outData=np.zeros((2,self.length))
        self.derivative=np.zeros(self.length)
        self.derivativeN=np.zeros(self.length)
        self.CriticalPoints=[]
        self.Peaks=[]
        self.Troughs=[]
        self.Inflections=[]
        self.IndividualWF=[]
        self.IndividualWFaverage=[]
        self.fs=60 #sampling rate, later determine from the data
        InitTime=self.ISO2mSec(self.rawData.pop(0))
        i=self.length-1;
        while self.rawData:
            oneData=self.rawData.pop()            
            self.Time[i]=self.ISO2mSec(oneData)-InitTime
            self.Data[i]=oneData[2]
            i=i-1

    def __init__(self, openF):
        self.make(openF)

    def getMaxTime(self):
        max=self.Time[0]
        for k in range(0,self.length):
            if self.Time[k]>min:
                max=self.Time
        return max
    
    def getWF(self):
        self.outData[0]=self.Time
        self.outData[1]=self.Data
        return np.copy(self.outData)

    def getNormalized(self, temp):
        templen=len(temp)
        self.tempData=np.zeros(templen)
        max=temp[0]
        for k in range(0,templen):
            if temp[k]>max:
                max=temp[k]
        for k in range(0,templen):
            self.tempData[k]=temp[k]/max
        #print "Normalized", self.tempData
        return np.copy(self.tempData)


    def WFnormP(self):
        self.getWFsmooth()
        self.outData[1]= self.getNormalized(self.outData[1])
        return np.copy(self.outData)

    def getWFnormP(self):
        self.WFnormP()
        self.DataSN=np.copy(self.outData[1])
        return np.copy(self.outData)
  
        
    def getWFsmooth(self):
        self.outData[0]=self.Time
        self.outData[1][0]=self.Data[0]
        for k in range(1,self.length-1):
            self.outData[1][k]=(self.Data[k-1]+self.Data[k]+self.Data[k+1])/3
        self.outData[1][self.length-1]=self.Data[self.length-1]            
        return np.copy(self.outData)

#    def getSlope(self,t0,t1,y0,y1):
#        deltatime=t1-t0
#        slope=0
#        if deltatime==0:
#            slope=float('nan')
#        else:
#            slope= (y1-y0)/deltatime
#        #print 'Time= ',t1,'    slope= ', slope
#        return slope
            
    def getWFderivative(self):
        sign=1
        #print 'Smoothed Data',self.DataSN
        self.WFnormP()
        for k in range(0,self.length-1):
            self.derivative[k]=PPGfunc.getSlope(self.outData[0][k],self.outData[0][k+1],self.outData[1][k],self.outData[1][k+1])
        self.derivative[self.length-1]=PPGfunc.getSlope(self.outData[0][self.length-2],self.outData[0][self.length-1],self.outData[1][self.length-2],self.outData[1][self.length-1])
        #print "DerivativeNow", self.derivative
        self.derivativeN=self.getNormalized(self.derivative)
        for k in range(0,self.length-1):
            if self.derivativeN[k]==float('nan'):
                self.derivativeN[k]=sign
            else:
                sign=self.derivativeN[k]/self.derivativeN[k]
        if self.derivativeN[self.length-1]==float('nan'):
            self.derivativeN[k]=sign
        #for k in range(0,self.length-1):
                #if self.derivativeN[k]==float('nan'):
                    #print 'NAN= ', k
        self.outData[1]=self.derivativeN
        return np.copy(self.outData)
        
    def getWFCriticalP(self):
        self.getWFderivative()
        k=3
        while k < self.length-3:
            #print "k= ",k, "Time= ",self.Time[k],"    Slope= ", self.derivativeN[k], 'Slope= ', self.derivativeN[k+1],'Ratio= ', self.derivativeN[k]/self.derivativeN[k+1]
            #if self.derivativeN[k]==float('nan'):
                #print 'NAN'
            if (self.derivativeN[k]/self.derivativeN[k+1])<=0:
                #print 'Critical Points change sign= ', k
                self.CriticalPoints.append(k+1)
                if(self.derivativeN[k]>0 and self.derivativeN[k+1]<0):
                    self.Peaks.append(k+1)
                else:
                    self.Troughs.append(k+1)
                k=k+2
            else:
                count=0
#                while self.derivativeN[k+count]<0 and abs(self.derivativeN[k+count])<0.006:
#                    #print "Zero SLope= ",k
#                    count=count+1
#                if count>5:
#                    self.CriticalPoints.append(k)
#                    self.CriticalPoints.append(k+count)
#                    self.Inflections.append([k,k+count])
#                    #print 'Critical Points zero= ', k
#                    k=k+count
            k=k+1
            
        CP=np.zeros((2,len(self.CriticalPoints)))
        i=0
        for k in self.CriticalPoints:
            CP[0][i]=self.Time[k]
            CP[1][i]=self.DataSN[k]
            i=i+1
        return np.copy(CP)
    
    def AverageWF(self):
        waveF=[]
        peak=0
        trough=0
        foundPeak=False
        FoundinitTrough=False
        initTrough=0
        pos=0
        #print self.CriticalPoints
        k=0
        while(k < len(self.CriticalPoints)-1):
            k=k+1
            pos=self.CriticalPoints[k]
            if self.Time[pos]>5000 and self.Time[pos]<15000:
                if FoundinitTrough==False:
                    try:
                        self.Troughs.index(pos)
                        FoundinitTrough=True
                        initTrough=pos
                    except ValueError:
                       continue
                elif foundPeak==False:
                    #print 'pos= ',pos, ' Time= ', self.Time[pos]
                    try:
                        self.Peaks.index(pos)
                        foundPeak=True
                        peak=pos
                        #print 'Peak= ',peak, ' Time= ', self.Time[pos]
                    except ValueError:
                        initTrough=pos
                        continue
                else:
                    try:
                        self.Troughs.index(pos)
                        if self.Time[pos]-self.Time[ peak]<200:
                            #print 'Skip Trough pos= ',pos, ' Time= ', self.Time[pos]
                            k=k+1
                            #print 'Skip Trough pos= ',self.CriticalPoints[k], ' Time= ', self.Time[self.CriticalPoints[k]]
                            #print 'Skip Trough pos= ',self.CriticalPoints[k+1], ' Time= ', self.Time[self.CriticalPoints[k+1]]
                        else:
                            foundPeak=False
                            self.IndividualWF.append([initTrough,peak,pos])
                            initTrough=pos
                            #print 'Troughs= ',initTrough, ' Time= ', self.Time[pos]
                    except ValueError:
                       continue
        #print "Waves = ", self.IndividualWF
        numWave=0
        lenWave=0;
        MaxlenWave=0;
        timeWave=0
        MaxTimeWave=0
        timePeak=0
        
        for k in self.IndividualWF:
            #print "k waves= ",k
            numWave=numWave+1
            lenWave=k[2]-k[0]
            timeWave=self.Time[k[2]]-self.Time[k[0]]
            timePeak=timePeak+self.Time[k[1]]-self.Time[k[0]]
            if timeWave>MaxTimeWave:
                MaxtimeWave=timeWave
            #print 'Wave from ',self.Time[k[0]],'to ', self.Time[k[2]]
            if lenWave>MaxlenWave:
                MaxlenWave=lenWave
            deltatime=MaxtimeWave/MaxlenWave
            #print 'Number of Waves= ',  numWave,'Wave length= ',  lenWave, ' Max Wave length= ', MaxlenWave, ' Max Wave Time= ',  MaxtimeWave, ' delta Time= ',  deltatime,  'Peak Time ', self.Time[k[1]]-self.Time[k[0]], ' Average peak Time= ',  timePeak

        timePeak=timePeak/numWave
        #print 'Avererage Peak Time', timePeak
        MaxLen=MaxlenWave #*2
        WaveA=np.zeros(MaxLen)
        WaveT=np.zeros(MaxLen)
        Count=np.zeros(MaxLen)
        AverageData=np.zeros((3,MaxLen))
        #deltatime=MaxtimeWave/(MaxLen)
        PeakTimePoint=int(timePeak/deltatime)
        #print 'Peak Time=' ,timePeak, 'Peak Time Point=', PeakTimePoint
        
        i=0
        for k in self.IndividualWF:

            i=PeakTimePoint+1
            
            minValue=self.Data[k[0]]
            for j in range(k[0], k[2]):
                if self.Data[j]<minValue:
                    minValue=self.Data[j]

            AverageData[1][PeakTimePoint]=AverageData[1][PeakTimePoint]+self.Data[k[1]]-minValue
            Count[PeakTimePoint]=Count[PeakTimePoint]+1
            #print  AverageData[1][PeakTimePoint],  AverageData[1][PeakTimePoint]/numWave

            for j in range(k[1]+1, k[2]):
                i=round((self.Time[j]-self.Time[k[0]])/deltatime)
                #print 'Index= ',i
                if i==PeakTimePoint:
                    i=PeakTimePoint+1
                if i>=MaxLen:
                    i=MaxLen-1
                AverageData[1][i]=AverageData[1][i]+self.Data[j]-minValue
                Count[i]=Count[i]+1
                i=i+1
            #print  'Average Data=', AverageData[1]
                
          
        for k in self.IndividualWF:
            i=0
            for j in range(k[0], k[1]-1):
                i=round((self.Time[j]-self.Time[k[0]])/deltatime)
                #print 'Index= ',i
                if i==PeakTimePoint:
                    i=PeakTimePoint-1
                AverageData[1][i]=AverageData[1][i]-minValue+self.Data[j]
                Count[i]=Count[i]+1

                i=i+1
            #print  'Average Data=', AverageData[1]

        for i in range(0,MaxLen):
            AverageData[0][i]=deltatime*i
            AverageData[1][i]=AverageData[1][i]/Count[i]
        #print 'Peak Time= ',  PeakTimePoint, 'Max Len= ', MaxLen
        AverageData[2][0]=PeakTimePoint
        AverageData[2][1]=MaxLen-1
        #print 'Peak Time= ',  AverageData[2][0], 'Max Len= ', AverageData[2][1]

        #print 'Index= ',i, 'Time=' ,AverageData[0][i], 'Average=', AverageData[1][i]


        #print 'Average Data', AverageData


        return np.copy(AverageData)


    def Spectrum(self):
        n = len(self.Data) # length of the signal
        k = arange(n)
        T = n/self.fs
        frq = k/T # two sides frequency range
        frq = frq[range(n/2)] # one side frequency range

        Y = fft(self.Data)/n # fft computing and normalization
        #print Y
        Y = Y[range(n/2)]
        #print Y

        x1 = frq
        y1 = abs(Y)
        plt.scatter(x1, y1,label = "Frequency Distibution", color= "red", marker= "*", s=10)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Distibution [Arb. Units]')
        plt.title('PGP Frequency Domain')
        plt.xlim([0,100])
        plt.ylim([0,50])
        plt.show()
       
#        plot(frq,abs(Y),'r') # plotting the spectrum
#        xlabel('Freq (Hz)')
#        ylabel('|Y(freq)|')

#        Fs = 150.0  # sampling rate
#        Ts = 1.0/Fs # sampling interval
#        t = arange(0,1,Ts) # time vector
#        ff = 5   # frequency of the signal
#        y = sin(2*pi*ff*t)
#        subplot(2,1,1)
#        plot(t,y)
#        xlabel('Time')
#        ylabel('Amplitude')
#        subplot(2,1,2)
#        plotSpectrum(y,Fs)
#        show()

        
    def getWFbandPass(self, low,high):
        for k in range(0,self.length):
            self.outData[0][i]=self.time[i]
            self.outData[1][i]=self.data[i]
        self.outData[2][0]=self.data[i]            
        return np.copy(self.outData)
        

