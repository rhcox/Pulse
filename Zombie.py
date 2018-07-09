import board
import neopixel
import time
import adafruit_thermistor
from analogio import AnalogIn
from simpleio import map_range
import digitalio

def fileerror(text,e):
    print(text,' ',e)
    valuel=(255,0,0)
    if e.args[0] == 30:
        valuel=(0,0,255)
    if e.args[0] == 28:
        valuel=(0,255,0)
    pixels.fill((0, 0, 0))
    pixels[9] = valuel
    pixels.show()

    
def getvalue(line,n):
    string=''
    comma=0
    for c in line:
        if(comma<n):
            if c==',':
                comma=comma+1
        else:
            if c==',':
                break
            else:
                string=string+c
    return float(string)
    
def gettime(line):
    return getvalue(line,0)
    
def gettemp(line):
    return getvalue(line,1)
    
def getsm(line):
    return getvalue(line,1)

def getpulse(line):
    return getvalue(line,2)
    
def smooth(fp):
    count=0
    ftime=[0,0,0,0,0]
    fpulse=[0,0,0,0,0]
    try:
        with open("zombiesm.csv", "w+",0) as fp1:
            for line in fp:
                if count>0:
                    ftime[count-1]=gettime(line)
                    fpulse[count-1]=getpulse(line)
                    if(count>4):
                        total=0
                        for i in fpulse:
                            total=total+i
                        fp1.write('{0:f},'.format(ftime[2]))
                        fp1.write('{0:f}\n'.format(total/5))
                        for i in range(0,4):
                            ftime[i]=ftime[i+1]
                            fpulse[i]=fpulse[i+1]
                        count=5
                    else:
                        count=count+1
                else:
                    count=count+1
    except OSError as e:
        fileerror('Zombie Smooth',e)

def pulseRate():
    pulse=0
    try:
        with open("zombiesm.csv", "r",0) as fp:
            line=fp.readline()
            smmin=getsm(line)
            smmax=smmin
            for line in fp:
                val=getsm(line)
                if(val<smmin):
                    smmin=val
                if(val>smmax):
                    smmax=val
        with open("zombiesm.csv", "r",0) as fp:
            cutoff=smmax-smmin
            cutoff=smmin+int(float(cutoff)*0.25)
            pulse=0
            flag=True
            for line in fp:
                val=getsm(line)
                if val<cutoff:
                    if flag:
                        pulse=pulse+1
                    flag=False
                else:
                    flag=True
    except OSError as e:
        fileerror('Zombie Smooth Read',e)
    return pulse


def bodyTemp():
    try:
        with open("zombie.csv", "r",0) as fp:
            line=fp.readline()
            line=fp.readline()            
            smmin=gettemp(line)
            smmax=smmin
            sumTemp=smmin
            count=1
            for line in fp:
                val=gettemp(line)
                sumTemp=sumTemp+val
                if(val<smmin):
                    smmin=val
                if(val>smmax):
                    smmax=val
                count=count+1
            if (smmax+smmin)/2 < (sumTemp/count):
                return False
    except OSError as e:
        fileerror('Zombie Body Temperature',e)
    return True

pixels = neopixel.NeoPixel(board.NEOPIXEL, 10, brightness=.2)

thermistor = adafruit_thermistor.Thermistor(board.TEMPERATURE, 10000, 10000, 25, 3950)


analogin = AnalogIn(board.LIGHT)

WaitTime=5.0
DataTime=60.0
deltat=0.05
pixelTime=int(float(DataTime-WaitTime)/5.0+.5)
try:
    with open("zombie.csv", "w+",0) as fp:
        fp.write('Time, Temperature, Pulse Wave\n')
        fp.flush()
        pixels.fill((0, 100, 0))
        pixels.show()
        j=0
        timew=0
        while timew<WaitTime:
            time.sleep(deltat)
            j=j+1
            timew=j*deltat
        pixels.fill((0, 0, 0))
        pixels.show()
        
        pixels[1] = (0,255,0)
        pixels.show()
        pixelON=5
        
        pixelStart=timew
        timed=0
        j=0
        while timed<DataTime:
            temp = thermistor.temperature
            intensity=float(analogin.value)
            #it takes a long time to write
            fp.write('{0:f},'.format(timed))
            fp.write('{0:f},'.format(temp))
            fp.write('{0:f}\n'.format(intensity))
            fp.flush()
            time.sleep(deltat)
            j=j+1
            #so add that time adding like this is bad
            timed=timew+j*(deltat+0.40)
            if (timed-pixelStart)>pixelTime:
                if pixelON>9:
                    pixelON=9
                pixels[pixelON] = (0,0,50)
                pixels.show()
                pixelON=pixelON+1
                pixelStart=timed
        pixels.fill((0, 0, 0))
        pixels.show()
        fp.flush()
except OSError as e:
    fileerror('Zombie Write',e)
    

pixels[4] = ((20, 0, 0))
pixels.show()

try:
    with open("zombie.csv", "r",0) as fp:
        smooth(fp)
except OSError as e:
    fileerror('Zombie Read',e)

pixels[3] = ((20, 0, 0))
pixels.show()

pulseR=pulseRate()
pixels[2] = ((20, 0, 0))
pixels.show()

pulseR=int(float(pulseR)/(float(DataTime-WaitTime)/60.0))

pixels.fill((0, 0, 0))
pixels.show()
pixelc=pulseR/3
if pixelc>9:
    pixelc=9
    
if bodyTemp():
    pixelc=3
    
for i in range(0,pixelc):
    if i<5:
        pixels[i] = ((0,0,100))
    elif i<8:
        pixels[i] = ((0,100,0))
    else:
        pixels[i] = ((100,0,0))

