import neopixel
import board
pixels = neopixel.NeoPixel(board.NEOPIXEL, 10, brightness=.2)
pixels[0] = ((255,255,255))
pixels.show()

pixels[1] = ((255,255,255))
pixels.show()

import digitalio
switch = digitalio.DigitalInOut(board.SLIDE_SWITCH)
switch.direction = digitalio.Direction.INPUT
switch.pull = digitalio.Pull.UP

if switch.value:
    print(1)
    color=(0,0,100)
else:
    print(0)
    color=(0,100,0)
pixels[2] = color
pixels.show()

import time
pixels[3] = color
pixels.show()

import storage

pixels[4] = color
pixels.show()
    
for i in range(5,9):
    pixels[i] = color
    pixels.show()
    time.sleep(0.05)

time.sleep(1)

pixels.fill((0, 0, 0))
pixels.show()

storage.remount("/", switch.value)
