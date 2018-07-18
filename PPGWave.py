import time
import math
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
from scipy import signal
#from scipy.signal import butter, sosfilt, sosfreqz

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
        self.ReWave=[]
        self.InterpWave=[]
        self.IndividualWFaverage=[]
        self.AverageWave=np.zeros((3,self.length))
        self.AverageWF()

        self.troughThesh=0
        self.peakThesh=0
        self.StdDevPWD=0
        self.StdDevPWA=0

        

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
        return self.maxValue

    def getMinValue(self):
        self.minValue=PPGfunc.getMin(self.Data)
        return self.minValue

    def getMaxFiltered(self):
        self.maxFiltered=PPGfunc.getMax(self.filtered_sig)
        return self.maxFiltered

    def getMinFiltered(self):
        self.minFiltered=PPGfunc.getMin(self.filtered_sig)
        return self.minFiltered

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
#        return np.copy(outData)


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
        self.filtered_fft[np.abs(frq) < 0.5] = 0
        self.filtered_sig = abs(fftpack.ifft(self.filtered_fft))

        return self.filtered_sig


    def fileteredWave(self):
        n = len(self.Data)/100
        if n<5:
            n=5
        n=5
        #print ('n=',n)
        lowcut=0.5
        #print ('lowcut=',lowcut)
        highcut=5
        #print ('highcut=',highcut)
        fs=60
        #print ('fs=',fs)
        b,a=PPGfunc.butter_bandpass(lowcut, highcut, fs, n)
        #print ('Butter')
        self.filtered_sig=PPGfunc.butter_bandpass_filter(self.Data,b,a)
        #print ('Filter')
        #self.PlotPPGWave()
        #self.PlotPPGWaveFiltered()
        min=self.filtered_sig[0]
        for i in self.filtered_sig:
            if i<min:
                min=i
        self.filtered_sig=self.filtered_sig-min
        #print self.filtered_sig[10]
        #self.PlotPPGWaveBoth()

            
    def getWFderivative(self):
        self.Spectrum()
        self.fileteredWave()
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
#            print 'Is Trough=',self.Time[pos],self.filtered_sig[pos],self.troughThesh
            if self.filtered_sig[pos]<self.troughThesh:
 #               print 'YES'
                return True
        #    if (self.filtered_sig[pos]-self.minFiltered)<0.5*(self.maxFiltered-self.minFiltered):
        #        if peak>0:
        #            if (self.Time[pos]-self.Time[ peak])>200:
        #                return True
        #        return True
        return False
        
    def isPeak(self,pos,name,peak,trough):
        if name=="Peak":
#            print 'Is Peak=',self.Time[pos],self.filtered_sig[pos],self.peakThesh
            if self.filtered_sig[pos]>self.peakThesh:
 #               print 'Yes'
                #print 'Signal= ',self.filtered_sig[pos], 'Peak Threshhold=',self.peakThesh
                return True
#            if (self.filtered_sig[pos]-self.minFiltered)>0.60*(self.maxFiltered-self.minFiltered):
#                if peak>0:
#                    if (self.Time[pos]-self.Time[ peak])>200:
#                        return True
#                else:
#                    return True
        return False
        
    def findPeaksTroughs(self):
        k=0
        peak=-1
        trough=-1
        foundPeak=False
        FoundinitTrough=False
        initTrough=-1
        count=0
        Notch=0
        #self.PlotPPGWave()
        self.troughThesh=self.minFiltered+(self.maxFiltered-self.minFiltered)*.3
        self.peakThesh=self.maxFiltered-(self.maxFiltered-self.minFiltered)*.3
        for k in self.CriticalPoints:
            pos=k[0]
            name=k[1]
            if self.Time[pos]>10000:
                count=count+1
                if count>7:
                    peak=-1
                    trough=-1
                    foundPeak=False
                    FoundinitTrough=False
                    initTrough=-1
                    count=0
                if FoundinitTrough==False:
                    #print 'looking for init trough'
                    if self.isTrough(pos,name,peak,trough):
                        FoundinitTrough=True
                        initTrough=pos
                        peak=-1
                        count=0
                        #print "Found Init Trough"
                elif foundPeak==False:
                    #print "Looking for Peak"
                    if self.isPeak(pos,name,peak,trough):
                        foundPeak=True
                        peak=pos
                        count=0
                    elif self.isTrough(pos,name,peak,trough):
                       initTrough=pos
                       peak=-1
                       count=0
                else:
                    if self.isTrough(pos,name,peak,trough):
                        foundPeak=False
                        #print 'Found Wave Form= ',initTrough,peak,pos
                        minv=self.filtered_sig[initTrough]
                        for i in range(initTrough, pos):
                            if self.filtered_sig[i]<minv:
                                minv=self.filtered_sig[i]
                        self.IndividualWF.append([initTrough,peak,pos,count,minv])
                        #print 'Add Waveform= ',len(self.IndividualWF),self.IndividualWF
                        initTrough=pos
                        peak=-1;
                        count=0

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
        removeWF1=[]
        lenWave.sort()
        lenWaveLow=lenWave[len(lenWave)/5]
        lenWaveHi=lenWave[4*len(lenWave)/5]
        timeWave.sort()
        timeWaveLow=timeWave[len(timeWave)/5]
        timeWaveHi=timeWave[4*len(timeWave)/5]
        AampWave=ampWave/numWave
        
        i=0
        for k in (self.IndividualWF):
            lenWave=k[2]-k[0]
            timeWave=self.Time[k[2]]-self.Time[k[0]]
            ampWave=max(self.Data[k[1]]-self.Data[k[0]],self.Data[k[1]]-self.Data[k[2]])
            if lenWave<lenWaveLow or lenWave>lenWaveHi:
                removeWF.append(k)
            elif timeWave<timeWaveLow or timeWave>timeWaveHi:
                removeWF.append(k)
            elif (k[1]-k[0])>math.ceil(lenWave*0.5):
                removeWF.append(k)
            i=i+1
        wavesN=0
        
        for k in removeWF:
            #print 'remove=',k
            self.IndividualWF.remove(k)


    
        for k in (self.IndividualWF):
            if k[3]>4:
                #print 'No Notch=',k
                wavesN=wavesN+1
        if wavesN>5:
            for k in (self.IndividualWF):
                if k[3]<5:
                    removeWF1.append(k)
        
        for k in removeWF1:
            self.IndividualWF.remove(k)

        times=np.zeros(len(self.IndividualWF))
        i=0
        for k in self.IndividualWF:
            times[i]=self.Time[k[2]]-self.Time[k[0]]
            i=i+1
        self.StdDevPWD=np.std(times)
        print 'Standard Deviation PWD=', i, self.StdDevPWD
        
        PWApeak=np.zeros(len(self.IndividualWF))
        i=0
        for k in self.IndividualWF:
            PWApeak[i]=self.filtered_sig[k[1]]-k[4]
            i=i+1
        self.StdDevPWA=np.std(PWApeak)
        print 'Standard Deviation PWA=', i, self.StdDevPWA

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
        sumTime=0
        for k in self.IndividualWF:
            numWave=numWave+1
            lenWave=k[2]-k[0]
            timeWave=self.Time[k[2]]-self.Time[k[0]]
            sumTime=sumTime+timeWave
            if timeWave>MaxTimeWave:
                MaxTimeWave=timeWave
            if lenWave>MaxlenWave:
                MaxlenWave=lenWave
            timePeak=timePeak+self.Time[k[1]]-self.Time[k[0]]
        MaxLen=MaxlenWave*2
        if MaxLen<10:
            MaxLen=10
            
        if MaxlenWave>0 and numWave>0:
            wTime=sumTime/numWave
            deltatime=wTime/MaxLen


        Count=np.zeros(MaxLen)
        self.AverageWave.resize((3,MaxLen))
        
        for k in range(0,MaxLen):
            for i in range(1,3):
                self.AverageWave[i][k]=0
            self.AverageWave[0][k]=deltatime*k
        if numWave==0:
            self.AverageWave[2][2]=-1
            return
        timePeak=timePeak/numWave
        PeakTimePoint=int(timePeak/deltatime)

        tempTime1=np.zeros(PeakTimePoint)
        for i in range(0, PeakTimePoint):
            tempTime1[i]=deltatime*i

        tempTime2=np.zeros(MaxLen-PeakTimePoint)
        for i in range(PeakTimePoint,MaxLen):
            tempTime2[i-PeakTimePoint]=deltatime*i

        for k in self.IndividualWF:
            count=0
            
            i=0
            tempTime=np.zeros(k[2]-k[0])
            tempWave=np.zeros(k[2]-k[0])
            for j in range(k[0],k[2]):
                tempTime[i]=self.Time[j]-self.Time[k[0]]
                tempWave[i]=self.filtered_sig[j]-k[4]
                i=i+1
            self.ReWave.append(signal.resample(tempWave, MaxLen))
            
            
            tempT1=np.zeros(PeakTimePoint)
            tempWave1=np.zeros(PeakTimePoint)
            for j in range(PeakTimePoint):
                tempT1[j]=self.AverageWave[0][j]
                tempWave1[j]=self.ReWave[count][j]
            InterpWave1=np.interp(tempTime1, tempT1, tempWave1)
            #print InterpWave1

            tempT2=np.zeros(MaxLen-PeakTimePoint)
            tempWave2=np.zeros(MaxLen-PeakTimePoint)
            i=0
            for j in range(PeakTimePoint,MaxLen):
                tempT2[i]=self.AverageWave[0][j]
                tempWave2[i]=self.ReWave[count][j]
                i=i+1
            InterpWave2=np.interp(tempTime2, tempT2, tempWave2)

            tempWave3=np.zeros(MaxLen)
            for i in range(0,PeakTimePoint):
                tempWave3[i]=InterpWave1[i]
            for i in range(PeakTimePoint,MaxLen):
                tempWave3[i]=InterpWave2[i-PeakTimePoint]
                
            self.InterpWave.append(tempWave3)


        
            for i in range(0, MaxLen):
                self.AverageWave[1][i]=self.AverageWave[1][i]+tempWave3[i]
                Count[i]=Count[i]+1
               
                
            

            self.AverageWave[1][PeakTimePoint]=self.AverageWave[1][PeakTimePoint]+self.filtered_sig[k[1]]-k[4]
            Count[PeakTimePoint]=Count[PeakTimePoint]+1

        for i in range(0,MaxLen):
            self.AverageWave[1][i]=self.AverageWave[1][i]/Count[i]
        self.AverageWave[2][0]=PeakTimePoint
        self.AverageWave[2][1]=MaxLen-1
        self.AverageWave[2][2]=self.StdDevPWD
        self.AverageWave[2][3]=len(self.IndividualWF)
        self.AverageWave[2][4]=self.StdDevPWA

        i=0
        for k in self.IndividualWF:
            self.PlotPPGWaveFit(k,i,True)
            i=i+1
            
        return np.copy(self.AverageWave)

    

    def SetupPlotPPG(self,y1):
        x1 = self.Time
        plt.scatter(x1, y1,label = "Pulse Wave Data", color= "red", marker= "*", s=10)
        plt.plot(x1, y1,label = "Pulse Wave Form", color= "red")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='lower right')
        #plt.xlim(0,2000)
        plt.xlim([min(0,math.floor(self.minTime)),math.ceil(self.maxTime/100)*100])
        

    def PlotPPGWave(self):
        y1 = self.Data
        self.SetupPlotPPG(y1)
        plt.ylim([min(0,math.floor(self.minValue)),math.ceil(self.maxValue/10)*10])
        #print 'MinMax= ', self.minValue, self.maxValue
        plt.title('PGP Wave Raw')
        plt.show()

    def PlotPPGWaveBoth(self):
        x2 = self.Time
        y2 = self.filtered_sig
        plt.plot(x2, y2,label = "Pulse Wave Form filtered", color= "blue")
        y1=np.zeros(self.length)
        for i in range(0,self.length):
            y1[i] = self.Data[i]-self.minValue
        self.SetupPlotPPG(y1)
#        plt.ylim([min(0,math.floor(self.minValue)),math.ceil(self.maxValue/10)*10])
        #print 'MinMax= ', self.minValue, self.maxValue
        plt.title('PGP Wave Raw and Filtered')
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


    def PlotPPGWaveFit(self, k,count,flag):
        time=np.zeros(k[2]-k[0])
        wave=np.zeros(k[2]-k[0])
        i=0
        
        for j in range(k[0],k[2]):
            time[i]=self.Time[j]-self.Time[k[0]]
            wave[i]=self.filtered_sig[j]
            i=i+1

        x1 = time
        y1 = wave
        plt.scatter(x1, y1,label = "Filtered", color= "red", marker= "*", s=10)
        x3 = self.AverageWave[0]
        y3 = self.ReWave[count]
        plt.scatter(x3, y3,label = "Resampled", color= "blue", marker= "*", s=15)
        x2 = self.AverageWave[0]
        y2 =self.InterpWave[count]
        plt.plot(x2, y2,label = "Interpolated", color= "green")
        if flag:
            x4 = self.AverageWave[0]
            y4 = self.AverageWave[1]
            plt.plot(x4, y4,label = "Average Wave", color= "red")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.xlim(0,self.AverageWave[0][len(self.AverageWave[0])-1])
        plt.legend(loc='upper right')
        plt.title('PGP Sampling')
        plt.show()
