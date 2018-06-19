import numpy as np

def getSlope(t0,t1,y0,y1):
    deltatime=t1-t0
    slope=0
    if deltatime==0:
        slope=float('nan')
    else:
        slope= (y1-y0)/deltatime
        #print 'Time= ',t1,'    slope= ', slope
    return slope

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
            
