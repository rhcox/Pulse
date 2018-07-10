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
import PPGWave
import PPGfunc

    
class PPGSeg(object):
    def make(self,PPGData):
        self.ppg=PPGData
#       self.smoothS=smoothSpan
#       self.smooth= [[0 for x in range(columns)] for y in range(rows)]
#       self.smoothN=[

        self.PulseWaveSystolicPeakPos=-1
        self.PulseWaveEndPos=-1
        self.DicroticNotchPos=-1
        self.DiastolicPeakPos=-1
        self.HalfHeight=-1
        self.HalfTimeLow=-1
        self.HalfTimeHigh=-1                        
        self.PulseWaveAmplitude=-1
        self.SystolicPhase=-1
        self.DiastolicPhase=-1
        self.PulseWaveDuration=-1
        self.PulsePropogationTme=-1
        self.PulseWidthTime=-1
        self.RiseTime=-1
        self.PulseInflectionPointAmplitude=-1
        self.RelectionIndex=-1
        self.PulseArea=-1
        self.PulseArea1=-1
        self.PulseArea2=-1
        self.DicroticNotchTime=-1
        self.DicroticNotchValue=-1
        self.InflectionPointAreaRatio=-1
        self.HalfHeightHiPos=-1
        self.HalfHeightLowPos=-1
        
        self.AverageWaveForm=self.ppg.getAverageWF()
        self.PulseWaveSystolicPeakPos=int(self.AverageWaveForm[2][0])
        self.PulseWaveEndPos=int(self.AverageWaveForm[2][1])
        firstD=PPGfunc.derivative(self.AverageWaveForm[0],self.AverageWaveForm[1])
        secondD=PPGfunc.derivative(self.AverageWaveForm[0],firstD)



        x1 = self.AverageWaveForm[0]
        y1 = self.AverageWaveForm[1]
        plt.plot(x1, y1,label = "pulse", color= "red")
        x2 = self.AverageWaveForm[0]
        y2 = firstD*10
        plt.plot(x2, y2,label = "First Derivative", color= "blue")
        x2 = self.AverageWaveForm[0]
        y2 = secondD*100
        plt.plot(x2, y2,label = "Second Derivative", color= "green")
        plt.xlabel('Time [ms]')
        plt.ylabel('Pulse Waveform [Arb. Units]')
        plt.legend(loc='upper right')
        plt.title('Pulse Derivatives')
        plt.show()


        self.HalfHeight=0.5*self.AverageWaveForm[1][self.PulseWaveSystolicPeakPos]
        #print 'Half Height= ',self.HalfHeight
        for i in range(0,self.PulseWaveSystolicPeakPos-1):
            if self.AverageWaveForm[1][i]>self.HalfHeight:
                DeltaX=self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]
                DeltaY=self.AverageWaveForm[1][i]-self.AverageWaveForm[1][i-1]
                DeltaV=self.HalfHeight-self.AverageWaveForm[1][i-1]
                self.HalfTimeLow=self.AverageWaveForm[0][i-1]+DeltaX*DeltaV/DeltaY
                self.HalfHeightLowPos=i
                break
        for i in range(self.PulseWaveSystolicPeakPos+1,self.PulseWaveEndPos-1):
            if self.AverageWaveForm[1][i]<self.HalfHeight:
                DeltaX=self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]
                DeltaY=-self.AverageWaveForm[1][i]+self.AverageWaveForm[1][i-1]
                DeltaV=self.HalfHeight-self.AverageWaveForm[1][i]
                self.HalfTimeHigh=self.AverageWaveForm[0][i-1]+DeltaX*DeltaV/DeltaY
                self.HalfHeightHiPos=i
                break



        Inflections=PPGfunc.FindNotch(self.AverageWaveForm[0],self.AverageWaveForm[1],self.PulseWaveSystolicPeakPos,self.HalfHeightHiPos,firstD,secondD)
        if(len(Inflections)>1):
            self.DiastolicPeakPos=Inflections.pop()
            self.DicroticNotchPos=Inflections.pop()
            self.DicroticNotchTime=self.AverageWaveForm[0][self.DicroticNotchPos]
            self.DicroticNotchValue=self.AverageWaveForm[1][self.DicroticNotchPos]

        else:
            self.DicroticNotchPos=-1
            self.DiastolicPeakPos=-1

                        
        self.PulseWaveAmplitude=self.AverageWaveForm[1][self.PulseWaveSystolicPeakPos]
        self.SystolicPhase=self.RiseTime=self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.DiastolicPhase=self.AverageWaveForm[0][self.PulseWaveEndPos]-self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.PulseWaveDuration=self.AverageWaveForm[0][self.PulseWaveEndPos]
        self.PulsePropogationTme=self.AverageWaveForm[0][self.DiastolicPeakPos]-self.AverageWaveForm[0][self.PulseWaveSystolicPeakPos]
        self.PulseWidthTime=self.HalfTimeHigh-self.HalfTimeLow
        self.PulseInflectionPointAmplitude=self.AverageWaveForm[1][self.DiastolicPeakPos]
        self.RelectionIndex=self.PulseInflectionPointAmplitude/self.PulseWaveAmplitude*100

        self.PulseArea1=0
        for i in range(1,self.DicroticNotchPos):
            self.PulseArea1=self.PulseArea1+abs((self.AverageWaveForm[1][i]-self.AverageWaveForm[1][i-1])*(self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]))

        self.PulseArea2=0
        for i in range(self.DicroticNotchPos,self.PulseWaveEndPos):
            self.PulseArea2=self.PulseArea2+abs((self.AverageWaveForm[1][i]-self.AverageWaveForm[1][i-1])*(self.AverageWaveForm[0][i]-self.AverageWaveForm[0][i-1]))
        
        self.PulseArea=self.PulseArea1+self.PulseArea2
        self.InflectionPointAreaRatio=self.PulseArea2/self.PulseArea1

        
      
    def PlotPulse(self):
        x1 = self.AverageWaveForm[0]
        y1 = self.AverageWaveForm[1]
        plt.scatter(x1, y1,label = "PPG Wave Form", color= "red", marker= "*", s=10)
        plt.plot(x1, y1,label = "PPG Wave Form", color= "red")
        plt.xlabel('Time [ms]')
        plt.xlim([-int(self.PulseWaveDuration/50),self.PulseWaveDuration])
        plt.ylim([-int(self.PulseWaveAmplitude/5),self.PulseWaveAmplitude*1.25])
        
        if self.DicroticNotchPos>0:
            text='Dicrotic Notch: '+str(int(self.DicroticNotchTime))+' ms'
            plt.annotate(text,
            xy=(self.AverageWaveForm[0][self.DicroticNotchPos], self.AverageWaveForm[1][self.DicroticNotchPos]),
            xytext=(100, 50),textcoords='offset points', ha='right', va='bottom',
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
            
            plt.annotate('Diastolic Peak: '+str(int(self.AverageWaveForm[0][self.DiastolicPeakPos]))+' ms',
            xy=(self.AverageWaveForm[0][self.DiastolicPeakPos], self.AverageWaveForm[1][self.DiastolicPeakPos]),
            xytext=(150, 10), textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->',
            connectionstyle='arc3,rad=0'))

        text='PW= '+ str(int(self.PulseWidthTime))+' ms'
        start=self.HalfTimeLow+(self.HalfTimeHigh-self.HalfTimeLow)/2.0
        distance=(self.HalfTimeHigh-self.HalfTimeLow)/2.0
        plt.arrow(start, self.HalfHeight, -distance, 0, width=0.05,length_includes_head=True)
        plt.arrow(start, self.HalfHeight, distance, 0, width=0.05,length_includes_head=True)
        plt.text(start,self.HalfHeight , text, ha='center', va='top')
        
        text='PWA= '+ str(int(self.PulseWaveAmplitude))
        start=self.PulseWaveAmplitude/2.0
        distance=start
        plt.arrow(-2,start,  -2,-distance,width=0.05,length_includes_head=True)
        plt.arrow(-2,start,  -2, distance,width=0.05 ,length_includes_head=True)
        plt.text(start,self.HalfHeight , text, ha='right', va='center', rotation='vertical')

        self.RiseTime

        text='Systolic=RT= '+ str(int(self.RiseTime))+' ms'
        start=self.RiseTime/2
        distance=start
        plt.arrow(start,-2,  -distance,0,width=0.05,length_includes_head=True)
        plt.arrow(start,-2,  distance,0,width=0.05,length_includes_head=True)
        plt.text(start,-2 , text, ha='center', va='top')

        
        text='Diastolic= '+ str(int(self.DiastolicPhase))+' ms'
        start=self.RiseTime+self.DiastolicPhase/2.0
        distance=self.DiastolicPhase/2.0
        plt.arrow(start,-2,  -distance,0,width=0.05,length_includes_head=True)
        plt.arrow(start,-2,  distance,0,width=0.05,length_includes_head=True)
        plt.text(start,-2 , text, ha='center', va='top')
        
        text='PPT= '+ str(int(self.PulsePropogationTme))+' ms'
        start=self.RiseTime+self.PulsePropogationTme/2.0
        distance=self.PulsePropogationTme/2.0
        plt.arrow(start, self.PulseWaveAmplitude, -distance,0, width=0.05,length_includes_head=True)
        plt.arrow(start, self.PulseWaveAmplitude, distance, 0, width=0.05,length_includes_head=True)
        plt.text(start,self.PulseWaveAmplitude , text, ha='center', va='top')

        text='PIP= '+ str(int(self.PulseInflectionPointAmplitude))
        start=self.DicroticNotchValue/2.0
        distance=start
        plt.arrow(self.DicroticNotchTime, start, 0, -distance, width=0.05,length_includes_head=True)
        plt.arrow(self.DicroticNotchTime, start, 0, distance, width=0.05,length_includes_head=True)
 #       plt.text(start,self.PulseWaveAmplitude , text, ha='center', va='top')
        text='Pulse Area 1= '+ str(int(self.PulseArea1))
        plt.text(self.DicroticNotchTime/1.5,self.DicroticNotchValue/3.0 , text, ha='center', va='top')
        text='Pulse Area 2= '+ str(int(self.PulseArea2))
        plt.text(self.DicroticNotchTime/3.0+self.DicroticNotchTime,self.DicroticNotchValue/3.0 , text, ha='center', va='top')

        text='Reflection Index= '+ str(int(self.RelectionIndex))+'%'
        plt.text(self.DicroticNotchTime/2+self.DicroticNotchTime,self.PulseWaveAmplitude*3.0/4.0 , text, ha='center', va='top')

        text='IPA= '+str(round(self.InflectionPointAreaRatio,2))
        plt.text(self.DicroticNotchTime/2.0+self.DicroticNotchTime,self.PulseWaveAmplitude*7.0/8.0 , text, ha='center', va='top')


        plt.title('Annoted Averge Wave Form')

        
#        plt.annotate(text,
#        xy=(self.HalfTimeLow+(-self.HalfTimeLow+ self.HalfTimeHigh)/2,self.HalfHeight), xytext=(100, -20),
#        textcoords='offset points', ha='right', va='bottom',arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
 
        plt.show()



    def __init__(self,PPG):
        self.make(PPG)

    def getPWBegin(self):
        return self.PWF[self.getPWBegin]

    def getPWSystolicP(self):
        return self.PWF[self.PWPWSystolicP]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWDicroticNotch]

    def getPWDicroticNotch(self):
        return self.PWF[self.PWDicroticNotch]
    
    def getPWDiastolicPeak(self):
        return self.PWF[self.PWDiastolicPeak]
    
    def getPWEnd(self):
        return self.PWF[self.PWEnd]

    def getPWamplitude(self):
        return self.PWF[self.PWSystolicP]-self.PWF[self.PWBegin]
        
    def getRiseTime(self):
        return self.PWF[self.PWSystolicP]-self.PWF[self.PWBegin]

    def getPWDuration(self):
        return self.PWF[self.PWEnd]-self.PWF[self.PWBegin]

    def getPWPropTime(self):
        return self.PWF[self.PWDiastolicPeak]-self.PWF[self.PWSystolicP]

    def PrintAllParameters(self):
        print 'Pulse Wave Systolic Peak Position= ', self.PulseWaveSystolicPeakPos
        print 'Pulse Wave End Position= ', self.PulseWaveEndPos
        print 'Dicrotic Notch Position Position= ', self.DicroticNotchPos
        print 'Diastolic Peak Position= ',self.DiastolicPeakPos
        print  'Half Height= ', self.HalfHeight
        print  'Half Time Low= ', self.HalfTimeLow
        print  'HalfTimeHigh= ', self.HalfTimeHigh                      
        print  'Pulse Wave Amplitude= ', self.PulseWaveAmplitude
        print  'Systolic Phase= ', self.SystolicPhase, ' ms'
        print  'Rise Time= ', self.RiseTime, ' ms'
        print  'Diastolic Phase= ', self.DiastolicPhase, ' ms'
        print  'Pulse Wave Duration= ', self.PulseWaveDuration, ' ms'
        print  'Interbeat Interval= ', self.PulseWaveDuration, ' ms'
        print  'Pulse Propogation Tme= ', self.PulsePropogationTme, ' ms'
        print  'Pulse Width Time= ', self.PulseWidthTime, ' ms'
        print  'Pulse Inflection Point Amplitude= ', self.PulseInflectionPointAmplitude
        print 'Relection Index= ',self.RelectionIndex
        print 'Pulse Area= ', self.PulseArea
        print 'Pulse Area1= ',  self.PulseArea1
        print 'Pulse Area2= ', self.PulseArea2

#plotly.tools.set_credentials_file(username='rhcox', api_key='OWm9FUQ3FgSdgiuKs8Cq')
#p1=PPGWave.PPGWave()
#wf=PPGSeg(p1)
#p1.Spectrum()
#wf.PrintAllParameters()
