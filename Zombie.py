import board
import neopixel
import time
import adafruit_thermistor
from analogio import AnalogIn
from simpleio import map_range
import digitalio

pixels = neopixel.NeoPixel(board.NEOPIXEL, 10, brightness=.2)

thermistor = adafruit_thermistor.Thermistor(board.TEMPERATURE, 10000, 10000, 25, 3950)

pixels.fill((100,0,0))
pixels.show()

analogin = AnalogIn(board.LIGHT)

WaitTime=10.0
DataTime=60.0
deltat=0.05
pixelTime=int(float(DataTime-WaitTime)/5.0+.5)
try:
    with open("zombie.cvs", "w+",0) as fp:
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
        print('END Wait')
        
        pixels[1] = (0,255,0)
        pixels.show()
        pixelON=5
        
        pixelStart=timew
        timed=0
        j=0
        while timed<DataTime:
#            print('Collect Data= ',j,timed,DataTime)
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
            timed=timew+j*(deltat+0.15)
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
    print(e)
    delay = 0.5
    valuel=(255,0,0)
    if e.args[0] == 30:
        valuel=(0,0,255)
    if e.args[0] == 28:
        valuel=(0,255,0)
    pixels.fill((0, 0, 0))
    pixels[9] = valuel
    pixels.show()

pixels[5] = (50,0,0)
pixels.show()
