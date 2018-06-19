import PPGWave
import PPGSeg
import sys

if len(sys.argv) >1:
    openF=sys.argv[1]
else:
    openF='ppg.csv'
    
print 'File Name: ',openF
p1=PPGWave.PPGWave(openF)
wf=PPGSeg.PPGSeg(p1)
p1.Spectrum()
wf.PlotPPG()
wf.PlotWaveForm()
wf.PrintAllParameters()
