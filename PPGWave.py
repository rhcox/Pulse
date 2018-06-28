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
from scipy import fftpack

class PPGWave(object):
    def ISO2mSec(self, rData):
        dateString=rData[4]
        hour=int(dateString[11:13])
        minute=int(dateString[14:16])
        second=int(dateString[17:19])
        msecond=int(dateString[20:23])
        return (hour*60*60+minute*60+second)*1000+msecond

    def OpenFile(self,openF):
         with open(openF) as f:
            ppgData = csv.reader(f)
            self.rawData = [raw for raw in ppgData]    

    def getData(self):
        oneData=self.rawData.pop(0)
        if oneData[0]=='beep':
            i=0
            InitTime=self.ISO2mSec(self.rawData.pop(0))
            while self.rawData:
                oneData=self.rawData.pop(0)
                t=self.ISO2mSec(oneData)
                if float(oneData[2])>0 and float(t)>0:
                    self.Time[i]=t-InitTime
                    self.Data[i]=oneData[2]
                    i=i+1
            self.Time=np.resize(self.Time,i)
            self.Data=np.resize(self.Data,i)
            self.length=i
#            i=self.length-2;
#            InitTime=self.ISO2mSec(self.rawData.pop(0))
#            while self.rawData:
#                oneData=self.rawData.pop()            
#                self.Time[i]=self.ISO2mSec(oneData)-InitTime
#                self.Data[i]=oneData[2]
#                i=i-1
        elif oneData[0]=='amped_pulseWaveform':
            i=0
            while self.rawData:
                oneData=self.rawData.pop(0)
                if float(oneData[2])>0 and float(oneData[0])>0:
                    self.Time[i]=float(oneData[2])*1000
                    self.Data[i]=float(oneData[0])
                    i=i+1
            self.Time=np.resize(self.Time,i)
            self.Data=np.resize(self.Data,i)
            self.length=i
        else:
            print 'Input Format not Recognized'

    def make(self, openF):
        self.rawData=[]
        self.OpenFile(openF)

        self.length=len(self.rawData)-1
        self.Time=np.zeros(self.length)
        self.Data=np.zeros(self.length)
        self.getData()
        #print self.Data
        
        self.maxTime=0
        self.getMaxTime()
        self.minTime=0
        self.getMinTime()
        self.maxValue=0
        self.getMaxValue()
        self.minValue=0
        self.getMinValue()
        #print 'Maximum and Minimums= ',self.maxTime, self.minTime, self.maxValue, self.minValue

        self.derivative=np.zeros(self.length)
        self.dataFFT=np.zeros(self.length)
        self.filtered_fft = []
        self.filtered_sig = []
        self.radius=3
        self.fs=60 #sampling rate, later determine from the data
        self.getWFderivative()

        self.maxFiltered=0
        self.getMaxFiltered()
        self.minFiltered=0
        self.getMinFiltered()
        
        self.CriticalPoints=[]
        self.getWFCriticalPoints()

        self.Inflections=[]
        self.IndividualWF=[]
        self.IndividualWFaverage=[]
        self.AverageWave=np.zeros((3,self.length))
        self.AverageWF()
        

        

    def __init__(self, openF):
        self.make(openF)

    def getAverageWF(self):
        return np.copy(self.AverageWave)

    def getMaxTime(self):
        self.maxTime=PPGfunc.getMax(self.Time)
        return self.maxTime

    def getMinTime(self):
        self.minTime=PPGfunc.getMin(self.Time)
        return self.minTime
    
    def getMaxValue(self):
        self.maxValue=PPGfunc.getMax(self.Data)
        return self.maxTime

    def getMinValue(self):
        self.minValue=PPGfunc.getMin(self.Data)
        return self.minTime

    def getMaxFiltered(self):
        self.maxFiltered=PPGfunc.getMax(self.filtered_sig)
        return self.maxTime

    def getMinFiltered(self):
        self.minFiltered=PPGfunc.getMin(self.filtered_sig)
        return self.minTime

    def getWF(self):
        outData[0]=self.Time
        outData[1]=self.Data
        return np.copy(outData)

#    def getNormalized(self, temp):
#        templen=len(temp)
#        self.tempData=np.zeros(templen)
#        max=temp[0]
#        for k in range(0,templen):
#            if temp[k]>max:
#                max=temp[k]
#        for k in range(0,templen):
#            self.tempData[k]=temp[k]/max
#       return np.copy(self.tempData)


#    def WFnormP(self):
#        self.getWFsmooth()
#        self.outData[1]= self.getNormalized(self.outData[1])
#        return np.copy(self.outData)

#    def getWFnormP(self):
#        self.WFnormP()
#        self.DataSN=np.copy(self.outData[1])
#        return np.copy(self.outData)
  
        
#    def getWFsmooth(self):
#        outData[0]=self.Time
#        outData[1]=abs(self.Spectrum())
#        self.outData[1][0]=self.Data[0]
#        for k in range(1,self.length-1):
#            self.outData[1][k]=(self.Data[k-1]+self.Data[k]+self.Data[k+1])/3
#        self.outData[1][self.length-1]=self.Data[self.length-1]            
        return np.copy(outData)


    def Spectrum(self):
        n = len(self.Data) # length of the signal
        k = arange(n)
        T = n/self.fs
        frq = k/T # two sides frequency range
        frq = frq[range(n/2)] # one side frequency range

        self.dataFFT = fft(self.Data) # fft computing and normalization
        #PPGfunc.PlotFTT(frq,self.dataFFT,n)
        pos_mask = np.where(frq > 0)
        freqs = frq[pos_mask]
        peak_freq = freqs[self.dataFFT[pos_mask].argmax()]
        
        self.filtered_fft = self.dataFFT.copy()
        self.filtered_fft[np.abs(frq) > 5] = 0
        self.filtered_fft[np.abs(frq) < 1] = 0
        self.filtered_sig = abs(fftpack.ifft(self.filtered_fft))

        return self.filtered_sig

            
    def getWFderivative(self):
        self.Spectrum()
        for k in range(0,self.length-1):
            self.derivative[k]=PPGfunc.getSlope(self.Time[k],self.Time[k+1],self.filtered_sig[k],self.filtered_sig[k+1])
        self.derivative[self.length-1]=PPGfunc.getSlope(self.Time[self.length-2],self.Time[self.length-1],self.filtered_sig[self.length-2],self.filtered_sig[self.length-1])
        sign=1
        for k in range(0,self.length-1):
            if self.derivative[k]==float('nan'):
                self.derivative[k]=sign
            else:
                sign=self.derivative[k]/self.derivative[k]
        if self.derivative[self.length-1]==float('nan'):
            self.derivative[k]=sign
        return np.copy(self.derivative)
        
    def getWFCriticalPoints(self):
        k=3
        while k < self.length-3:
            if (self.derivative[k]/self.derivative[k+1])<=0:
                self.CriticalPoints.append([k+1,'Peak'])
                if(self.derivative[k]>0 and self.derivative[k+1]<0):
                    self.CriticalPoints.append([k+1,'Peak'])
                else:
                    self.CriticalPoints.append([k+1,'Trough'])
                k=k+2
            else:
                count=0
            k=k+1            
        CP=np.zeros((2,len(self.CriticalPoints)))
        i=0
        for k in self.CriticalPoints:
            CP[0][i]=self.Time[k[0]]
            CP[1][i]=self.filtered_sig[k[0]]
            i=i+1
        return np.copy(CP)

    def PlotCP(self):
        x1 = self.Time
        y1 = self.filtered_sig
        plt.plot(x1, y1,label = "Pulse Wave Form", color= "red")
        i=0
        x2 = np.zeros(len(self.CriticalPoints))
        y2 = np.zeros(len(self.CriticalPoints))
        for k in self.CriticalPoints:
            #print self.Time[k[0]],'    ',self.filtered_sig[k[0]]
            x2[i] = self.Time[k[0]]
            y2[i] = self.filtered_sig[k[0]]
            i=i+1
        plt.scatter(x2, y2,label = "Critical Points", color= "green", marker= "*", s=20)
        x3 = self.Time
        y3 = self.derivative*10
        plt.plot(x3, y3,label = "Derivative", color= "blue")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='lower right')
        plt.title('PGP Critical Points')
        plt.xlim([self.minTime,self.maxTime])
        plt.ylim([-5, self.maxFiltered])
        plt.show()

    def isTrough(self,pos,name,peak,trough):
        if name=="Trough":
            if (self.filtered_sig[pos]-self.minFiltered)<0.5*(self.maxFiltered-self.minFiltered):
                if peak>0:
                    if (self.Time[pos]-self.Time[ peak])>200:
                        return True
                return True
        return False
        
    def isPeak(self,pos,name,peak,trough):
        if name=="Peak":
            if (self.filtered_sig[pos]-self.minFiltered)>0.60*(self.maxFiltered-self.minFiltered):
                if peak>0:
                    if (self.Time[pos]-self.Time[ peak])>200:
                        return True
                else:
                    return True
        return False
        
    def findPeaksTroughs(self):
        
        k=0
        peak=-1
        trough=-1
        foundPeak=False
        FoundinitTrough=False
        initTrough=-1
        #self.PlotPPGWave()
        for k in self.CriticalPoints:
            #print 'Critical Point=',k
            pos=k[0]
            name=k[1]
            if self.Time[pos]>10000:
                if FoundinitTrough==False:
                    #print 'looking for init trough'
                    if self.isTrough(pos,name,peak,trough):
                        FoundinitTrough=True
                        initTrough=pos
                        peak=-1
                        #print "Found Init Trough"
                elif foundPeak==False:
                    #print "Looking for Peak"
                    if self.isPeak(pos,name,peak,trough):
                        foundPeak=True
                        peak=pos
                    elif self.isTrough(pos,name,peak,trough):
                       initTrough=pos
                       peak=-1
                else:
                    if self.isTrough(pos,name,peak,trough):
                        foundPeak=False
                        #print 'Found Wave Form= ',initTrough,peak,pos
                        self.IndividualWF.append([initTrough,peak,pos])
                        #print 'Add Waveform= ',len(self.IndividualWF),self.IndividualWF
                        initTrough=pos
                        peak=-1;
                        
        numWave=0
        lenWave=[]
        lenTop=0
        timeWave=[]
        ampWave=0
        for k in self.IndividualWF:
            numWave=numWave+1
            lenWave.append(k[2]-k[0])
            timeWave.append(self.Time[k[2]]-self.Time[k[0]])
            ampWave=max(self.Data[k[1]]-self.Data[k[0]],self.Data[k[1]]-self.Data[k[2]])
            
        removeWF=[]
        lenWave.sort()
        lenWaveLow=lenWave[len(lenWave)/5]
        lenWaveHi=lenWave[4*len(lenWave)/5]
        timeWave.sort()
        timeWaveLow=timeWave[len(timeWave)/5]
        timeWaveHi=timeWave[4*len(timeWave)/5]
        AampWave=ampWave/numWave
        
        i=0
        #print 'Length= ',len(self.IndividualWF), 'Before: ',self.IndividualWF
        for k in (self.IndividualWF):
            lenWave=k[2]-k[0]
            timeWave=self.Time[k[2]]-self.Time[k[0]]
            ampWave=max(self.Data[k[1]]-self.Data[k[0]],self.Data[k[1]]-self.Data[k[2]])
#            #print 'i=',i,self.Time[k[0]],self.Time[k[2]],'Len Wave=',lenWaveLow,lenWaveHi,lenWave,' Time Wave=',timeWaveLow,timeWaveHi,timeWave, 'Top Wave=',lenWave,k[1]-k[0]
            if lenWave<lenWaveLow or lenWave>lenWaveHi:
                removeWF.append(k)
                #print 'i=',i,lenWave,lenWaveLow,lenWaveHi,'len Wave Pop k', k
            elif timeWave<timeWaveLow or timeWave>timeWaveHi:
                removeWF.append(k)
                #print 'i=',i,timeWave,timeWaveLow,timeWaveHi,'time Wave Pop k', k
            elif (k[1]-k[0])>math.ceil(lenWave*0.5):
                removeWF.append(k)
                #print 'i=',i,k[1]-k[0],lenWave*0.5,'Top Wave Pop k', k
       #     elif ampWave>(AampWave*1.25) or ampWave<(AampWave*0.75):
        #        removeWF.append(k)
         #       print 'i=',i,'Amp Wave Pop k', k
#            self.PlotOnePPGWaveSeg(k,False)
            i=i+1
        for k in removeWF:
            self.IndividualWF.remove(k)

        #print 'Length= ',len(self.IndividualWF), 'After: ',self.IndividualWF
        return 1
                   
    def AverageWF(self):            
        waveF=[]
        self.findPeaksTroughs()

        #print "Waves = ", self.IndividualWF
        numWave=0
        lenWave=0;
        MaxlenWave=0;
        timeWave=0
        MaxTimeWave=0
        timePeak=0
        
        for k in self.IndividualWF:
            numWave=numWave+1
            lenWave=k[2]-k[0]
            timeWave=self.Time[k[2]]-self.Time[k[0]]
            if timeWave>MaxTimeWave:
                MaxTimeWave=timeWave
            if lenWave>MaxlenWave:
                MaxlenWave=lenWave
            timePeak=timePeak+self.Time[k[1]]-self.Time[k[0]]
            #print 'From ',self.Time[k[0]], ' to ',self.Time[k[2]]
        if MaxlenWave!=0:
            deltatime=MaxTimeWave/MaxlenWave

        
        MaxLen=MaxlenWave
        WaveA=np.zeros(MaxLen)
        WaveT=np.zeros(MaxLen)
        Count=np.zeros(MaxLen)
        self.AverageWave.resize((3,MaxLen))
        
        for k in range(0,MaxLen):
            for i in range(0,3):
                self.AverageWave[i][k]=0
        if numWave==0:
            self.AverageWave[2][2]=-1
            return

        timePeak=timePeak/numWave
        PeakTimePoint=int(timePeak/deltatime)
        #print 'Peak Time Point= ',PeakTimePoint
        

        for k in self.IndividualWF:
            i=PeakTimePoint+1
            minValue=self.filtered_sig[k[0]]
            for j in range(k[0], k[2]):
                if self.filtered_sig[j]<minValue:
                    minValue=self.filtered_sig[j]
            #print 'min Value', minValue

            self.AverageWave[1][PeakTimePoint]=self.AverageWave[1][PeakTimePoint]+self.filtered_sig[k[1]]-minValue
            Count[PeakTimePoint]=Count[PeakTimePoint]+1
            #print 'Accumulated Peak=',self.AverageWave[1][PeakTimePoint], 'Count=',Count[PeakTimePoint]

            for j in range(k[1]+1, k[2]):
                i=round((self.Time[j]-self.Time[k[0]])/deltatime)
                if i==PeakTimePoint:
                    i=PeakTimePoint+1
                if i>=MaxLen:
                    i=MaxLen-1
                self.AverageWave[1][i]=self.AverageWave[1][i]+self.filtered_sig[j]-minValue
                Count[i]=Count[i]+1
                i=i+1
                
          
        for k in self.IndividualWF:
            i=0
            for j in range(k[0], k[1]-1):
                i=round((self.Time[j]-self.Time[k[0]])/deltatime)
                #print 'Index= ',i
                if i==PeakTimePoint:
                    i=PeakTimePoint-1
                if i>=MaxLen:
                    i=MaxLen-1
                self.AverageWave[1][i]=self.AverageWave[1][i]-minValue+self.filtered_sig[j]
                Count[i]=Count[i]+1

                i=i+1

        for i in range(0,MaxLen):
            self.AverageWave[0][i]=deltatime*i
            self.AverageWave[1][i]=self.AverageWave[1][i]/Count[i]
        self.AverageWave[2][0]=PeakTimePoint
        self.AverageWave[2][1]=MaxLen-1

        return np.copy(self.AverageWave)

    

    def SetupPlotPPG(self,y1):
        x1 = self.Time
        plt.scatter(x1, y1,label = "Pulse Wave Data", color= "red", marker= "*", s=10)
        plt.plot(x1, y1,label = "Pulse Wave Form", color= "red")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='lower right')
        plt.xlim([min(0,math.floor(self.minTime)),math.ceil(self.maxTime/100)*100])
        

    def PlotPPGWave(self):
        y1 = self.Data
        self.SetupPlotPPG(y1)
        plt.ylim([min(0,math.floor(self.minValue)),math.ceil(self.maxValue/10)*10])
        #print 'MinMax= ', self.minValue, self.maxValue
        plt.title('PGP Wave Raw')
        plt.show()

    def PlotPPGWaveFiltered(self):
        y1 = self.filtered_sig
        self.SetupPlotPPG(y1)
        plt.title('PGP Wave Filtered')
        plt.show()

        
    def PlotOnePPGWaveSeg(self,k,flag):
        i=0
 
        if flag:
            x2=np.zeros(len(self.AverageWave[1]))
            y2=np.zeros(len(self.AverageWave[1]))
            y2=self.AverageWave[1]
            dtime=(self.Time[k[2]]-self.Time[k[0]])/float(len(x2))
            for i in range(0,len(self.AverageWave[1])):
                x2[i]=self.Time[k[0]]+i*dtime
            plt.plot(x2, y2,label = "Average Wave Form", color= "blue")


        x1=np.zeros(k[2]-k[0])
        y1=np.zeros(k[2]-k[0])
        j=0
        for i in range(k[0],k[2]):                               
            x1[j]=self.Time[i]
            y1[j]=self.filtered_sig[i]
            j=j+1
        plt.scatter(x1, y1,label = "Pulse Wave Form", color= "red", marker= "*", s=10)
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        time1=math.floor(math.floor(self.Time[k[0]]))
        time2=math.ceil(self.Time[k[2]]/10)*10
        plt.xlim([time1,time2])
        Text='PGP Wave Segmented from '+str(self.Time[k[0]])+' ms tto '+str(self.Time[k[2]])+'ms'
        plt.title(Text)
        plt.show()


    def PlotPPGWaveSeg(self):
#        i=0
#        x2=np.zeros(len(self.AverageWave[1]))
#        y2=np.zeros(len(self.AverageWave[1]))
#        y2=self.AverageWave[1]
        for k in self.IndividualWF:
            self.PlotOnePPGWaveSeg(k,True)
     #       dtime=(self.Time[k[2]]-self.Time[k[0]])/float(len(x2))
      #      #print 'dtime=',dtime
       #     for i in range(0,len(self.AverageWave[1])):
#                x2[i]=self.Time[k[0]]+i*dtime
#                #print 'i=',i,' dtime=',dtime,' i*dtime=',i*dtime,' x2[i]=', x2[i]
 #          #print 'X2= ',x2
  #          #print 'Y2= ',y2
#            plt.plot(x2, y2,label = "Average Wave Form", color= "blue")
  #          x1=np.zeros(k[2]-k[0])
 #           y1=np.zeros(k[2]-k[0])
  #          j=0
   #         for i in range(k[0],k[2]):                               
    #            x1[j]=self.Time[i]
 #               y1[j]=self.filtered_sig[i]
  #              j=j+1
   #         plt.scatter(x1, y1,label = "Pulse Wave Form", color= "red", marker= "*", s=10)
    #        plt.xlabel('Time [ms]')
 #           plt.ylabel('Pulse Waveform [Arb. Units]')
  #          time1=math.floor(math.floor(self.Time[k[0]]))
   #         time2=math.ceil(self.Time[k[2]]/10)*10
    #        #print 'Time 1= ',time1,'Time 2= ',time2
  #          plt.xlim([time1,time2])
   #         Text='PGP Wave Segmented from '+str(self.Time[k[0]])+' ms tto '+str(self.Time[k[2]])+'ms'
 #           plt.title(Text)
  #          plt.show()
