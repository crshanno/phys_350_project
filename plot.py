import os
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.realpath(__file__))
DT = 100 # seconds/frame

def importFile(filename):
    hndl = open(os.path.join(here, filename))

    time = []
    nums = []
    sizes = []
    sizeAvgs = []

    def processLine(line):
        if len(line) < 1:
            return
        line = line.strip()
        if line[0] == '#':
            return
        data = line.split(",")
        if len(data) < 2:
            print("Invalid line format")
            return
        t = int(data[0]) * DT
        n = int(data[1])
        del data[0]
        del data[1]
        for i in range(len(data)):
            data[i] = float(data[i])
        time.append(t)
        nums.append(n)
        sizes.append(data)
        sizeAvgs.append(sum(data)/len(data))

    buffersize = 2**16
    while True:
        lines = hndl.readlines(buffersize)
        if not lines:
            break
        for line in lines:
            processLine(line)
    hndl.close()

    plt.subplot(1,2,1)
    plt.plot(time, nums)
    plt.ylabel("Number of Objects")
    plt.xlabel("Time (s)")
    plt.subplot(1,2,2)
    plt.plot(time, sizeAvgs)
    plt.ylabel("Avg Size (m)")
    plt.xlabel("Time (s)")
    plt.show()
    
FILE = "1.csv"
importFile(FILE)
    
    