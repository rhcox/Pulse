import PPGWave
import PPGSeg
import sys

if len(sys.argv) >1:
    openF=sys.argv[1]
else:
    openF='ppg.csv'
    
print 'File Name: ',openF
p1=PPGWave.PPGWave(openF)
#p1.PlotPPGWave()
#p1.PlotPPGWaveBoth()
#p1.PlotPPGWaveFiltered()
#p1.PlotCP()
#p1.PlotPPGWaveSeg()

wf=PPGSeg.PPGSeg(p1)
wf.PlotPulse()
wf.PrintAllParameters()
