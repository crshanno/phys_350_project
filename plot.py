import os
import sys
import matplotlib.pyplot as plt
import math
import numpy

here = os.path.dirname(os.path.realpath(__file__))

def importFile(filename,label):
    path = (os.path.join(here, filename))

    data = numpy.genfromtxt(path, delimiter=',', skip_header=6, usecols=(0,1))
    
    time = data[:,0]
    nums = data[:,1]
    #totalMass = []

    #def processLine(line):
    #    if len(line) < 1:
    #        return
    #    line = line.strip()
    #    if line[0] == '#':
    #        return
    #    data = line.split(",")
    #    if len(data) < 2:
    #        print("Invalid line format")
    #        return
    #    t = float(data[0])
    #    n = int(data[1])
    #    del data[0]
    #    del data[1]
    #    for i in range(len(data)):
    #        data[i] = float(data[i])
    #    time.append(t)
    #    nums.append(n)
    #    sizes.append(data)
    #    sizeAvgs.append(sum(data)/len(data))
    #    
    #    mases = []
    #    for d in data:
            #d = 4./3. * math.pi * math.pow(d,3)
            #d = d * 132;
    #        mases.append(d)
    #    totalMass.append(sum(mases))

    #buffersize = 2**16
    #while False:
    #    lines = hndl.readlines(buffersize)
    #    if not lines:
    #        break
    #    for line in lines:
    #        processLine(line)
    #hndl.close()

    plt.subplot(121)
    h, = plt.loglog(time, nums, label=label)
    plt.subplot(122)
    plt.plot(time, nums, label=label)
    #plt.subplot(1,3,2)
    #plt.plot(time, sizeAvgs, label=label)
    #plt.ylabel("Average Size (m)")
    #plt.xlabel("Time (h)")
    #plt.subplot(1,3,3)
    #plt.plot(time, totalMass, label=label)
    #plt.ylabel("Total Mass (kg)")
    #plt.xlabel("Time (h)")
    print("done  %s" %filename)
    sys.stdout.flush()
    return h
    

h = []
def add(fn, lbl):
    h.append(importFile(fn,lbl))

if len(sys.argv) > 1:
    for i in range(1, len(sys.argv)):
        label = sys.argv[i]
        while label[0] == "_":
            label = label[1:]
        add(sys.argv[i], label)
else:
    pass

labelsize = 12
plt.subplot(121)
plt.grid(True)
plt.ylabel("Number of Objects", fontsize=labelsize)
plt.xlabel("Time (s)", fontsize=labelsize)
plt.title("")
plt.legend(handles=h)
axes = plt.gca()
#axes.set_xlim([0,10000])
plt.subplot(122)
plt.grid(True)
plt.ylabel("Number of Objects", fontsize=labelsize)
plt.xlabel("Time (s)", fontsize=labelsize)
plt.title("")
plt.legend(handles=h)
axes = plt.gca()
#axes.set_xlim([0,10000])
plt.show()

    