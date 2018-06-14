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
        
    def make(self):
        with open("ppg.csv") as f:
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

        InitTime=self.ISO2mSec(self.rawData.pop(0))
        i=self.length-1;
        while self.rawData:
            oneData=self.rawData.pop()            
            self.Time[i]=self.ISO2mSec(oneData)-InitTime
            self.Data[i]=oneData[2]
            i=i-1

    def __init__(self):
        self.make()

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

    def getSlope(self,t0,t1,y0,y1):
        deltatime=t1-t0
        slope=0
        if deltatime==0:
            slope=float('nan')
        else:
            slope= (y1-y0)/deltatime
        #print 'Time= ',t1,'    slope= ', slope
        return slope
            
    def getWFderivative(self):
        sign=1
        #print 'Smoothed Data',self.DataSN
        self.WFnormP()
        for k in range(0,self.length-1):
            self.derivative[k]=self.getSlope(self.outData[0][k],self.outData[0][k+1],self.outData[1][k],self.outData[1][k+1])
        self.derivative[self.length-1]=self.getSlope(self.outData[0][self.length-2],self.outData[0][self.length-1],self.outData[1][self.length-2],self.outData[1][self.length-1])
        #print "DerivativeNow", self.derivative
        self.derivativeN=self.getNormalized(self.derivative)
        for k in range(0,self.length-1):
            if self.derivativeN[k]==float('nan'):
                self.derivativeN[k]=sign
            else:
                sign=self.derivativeN[k]/self.derivativeN[k]
        if self.derivativeN[self.length-1]==float('nan'):
            self.derivativeN[k]=sign
        for k in range(0,self.length-1):
            if self.derivativeN[k]==float('nan'):
                print 'NAN= ', k
        self.outData[1]=self.derivativeN
        return np.copy(self.outData)
        
    def getWFCriticalP(self):
        self.getWFderivative()
        k=3
        while k < self.length-3:
            print "k= ",k, "Time= ",self.Time[k],"    Slope= ", self.derivativeN[k], 'Slope= ', self.derivativeN[k+1],'Ratio= ', self.derivativeN[k]/self.derivativeN[k+1]
            if self.derivativeN[k]==float('nan'):
                print 'NAN'
            if (self.derivativeN[k]/self.derivativeN[k+1])<=0:
                print 'Critical Points change sign= ', k
                self.CriticalPoints.append(k+1)
                if(self.derivativeN[k]>0 and self.derivativeN[k+1]<0):
                    self.Peaks.append(k+1)
                else:
                    self.Troughs.append(k+1)
                k=k+2
            else:
                count=0
                while self.derivativeN[k+count]<0 and abs(self.derivativeN[k+count])<0.006:
                    print "Zero SLope= ",k
                    count=count+1
                if count>5:
                    self.CriticalPoints.append(k)
                    self.CriticalPoints.append(k+count)
                    self.Inflections.append([k,k+count])
                    #print 'Critical Points zero= ', k
                    k=k+count
            k=k+1
            
        CP=np.zeros((2,len(self.CriticalPoints)))
        i=0
        for k in self.CriticalPoints:
            CP[0][i]=self.Time[k]
            CP[1][i]=self.DataSN[k]
            i=i+1
        return np.copy(CP)

        
    def getWFbandPass(self, low,high):
        for k in range(0,self.length):
            self.outData[0][i]=self.time[i]
            self.outData[1][i]=self.data[i]
        return np.copy(self.outData)
        



    
class PulseWave(object):
    def make(self,PPGData):
        self.ppg=PPGData
#       self.smoothS=smoothSpan
#       self.smooth= [[0 for x in range(columns)] for y in range(rows)]
#       self.smoothN=[]
        self.PWF=[[1.0,2.0,3.0,4.0,5.0], [0.0,0.2,0.4,0.6,0.8], [10.0,11.0,12.0,13.0,14.0]]
        self.wf1=self.ppg.getWFnormP()
        #print "Waveform0", self.wf1[1]
        #data_obj = Data(self.wf1[0], self.wf1[1], smoothness=11)
        #print "Waveform1", self.wf1[1]
        #data_obj.get_peaks(method='slope')
        #print "Waveform2", self.wf1[1]
        self.wfD=self.ppg.getWFderivative()
        self.wfCP=self.ppg.getWFCriticalP()
        #print "Waveform3", self.wf1[1]
        

        #print data_obj.peaks
        #print data_obj.peaks['peaks'][0]
        #print data_obj.peaks['peaks'][1]
        #print data_obj.peaks['valleys']

        #print 'Critical Points',self.wfCP


#        Pulse = go.Scatter(x=self.wf1[0],y=self.wf1[1],mode = 'markers',marker = dict(size = 1,color = 'rgba(255, 0, 0, .8)'),name='Pulse',)
#        Trace4 = go.Scatter(x=self.wfCP[0],y=self.wfCP[1],mode = 'markers',marker = dict(size = 3,color = 'rgba(255, 0, 0, .8)'),name='Pulse',)
        #Valleys = go.Scatter(x=data_obj.peaks['peaks'][0],y=data_obj.peaks['peaks'][1],mode = 'markers',marker = dict(size = 5,color = 'rgba(0, 255, 100, .8)'),name='Valleys',)
        #Peaks = go.Scatter(x=data_obj.peaks['valleys'][0],y=data_obj.peaks['valleys'][1],mode = 'markers',marker = dict(size = 5,color = 'rgba(0, 0, 255, .8)'),name='Peaks',)
#        data = [Pulse,Trace4]
#        layout = go.Layout(title='PPG',hovermode='closest',
#        xaxis=dict(title='Time [ms]',),
#        yaxis=dict(title='Pulse Waveform [Arb. Units]',),)
#        fig = go.Figure(data=data, layout=layout)
        #print data
#        py.plot(fig, filename='basic-scatter')
        #print "Out"

        x1 = self.wf1[0]
        y1 = self.wf1[1]
        plt.scatter(x1, y1,label = "Pulse Wave Form", color= "red", marker= "*", s=1)
        x2 = self.wfCP[0]
        y2 = self.wfCP[1]
        plt.scatter(x2, y2,label = "Critical Points", color= "green", marker= "*", s=20)
        x3 = self.wfD[0]
        y3 = self.wfD[1]
        plt.plot(x3, y3,label = "Derivative", color= "red")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='lower right')
        plt.title('PGP Critical Points')
        plt.xlim([7500,8500])
        plt.ylim([-.1,1])
        plt.show()

        
        #help(Data)
        self.PWTime = 0;
        self.PWNorm = 1;
        self.PWData = 2;       
        self.PWBegin = 0;
        self.PWSystolicP=1
        self.PWDicroticNotch=2
        self.PWDiastolicPeak=3
        self.PWEnd=4
        

    def __init__(self,PPG):
        self.make(PPG)

    def getPWBegin(self):
        return self.PWF[self.PWData][self.getPWBegin]

    def getPWSystolicP(self):
        return self.PWF[self.PWData][self.PWPWSystolicP]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWData][self.PWDicroticNotch]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWData][self.PWDicroticNotch]
    
    def getPWDiastolicPeak(self):
        return self.PWF[self.PWData][self.PWDiastolicPeak]
    
    def getPWEnd(self):
        return self.PWF[self.PWData][self.PWEnd]

    def getPWamplitude(self):
        return self.PWF[self.PWData][self.PWSystolicP]-self.PWF[self.PWData][self.PWBegin]
        
    def getRiseTime(self):
        return self.PWF[self.PWTime][self.PWSystolicP]-self.PWF[self.PWTime][self.PWBegin]

    def getPWDuration(self):
        return self.PWF[self.PWTime][self.PWEnd]-self.PWF[self.PWTime][self.PWBegin]

    def getPWPropTime(self):
        return self.PWF[self.PWTime][self.PWDiastolicPeak]-self.PWF[self.PWTime][self.PWSystolicP]

#plotly.tools.set_credentials_file(username='rhcox', api_key='OWm9FUQ3FgSdgiuKs8Cq')
p1=PPGWave()
wf=PulseWave(p1)
print wf.getPWamplitude()
print wf.getRiseTime()
print wf.getPWDuration()
print wf.getPWPropTime()
