import random
import math
import time
import os
import sys
import argparse
import numpy as np


#Use collision interpolation (check for collision between time steps)
useinterp = True

random.seed(time.time())

G = 6.67408*math.pow(10, -11) #m^3kg^-1s^-2 (Gravitational constant)
Mearth = 5.972*math.pow(10, 24) #kg (mass of earth)

Re = 6371 #km (radius of earth)

vstart=2
rstart = 1

pLow = 88*60
pHigh = 127*60

massMin = 1000 #kg
massMax = 10000 #kg

CollideScale = 10000 #scales likelyhood of collisions

AverageMass = 100
MinMass = 1
Density = 132 #kg/m^3 - based on https://en.wikipedia.org/wiki/Envisat

startBreakup = True
breakupTime = 1000
cb = []    
current_frame = 0
end_frame = 100

DT = 100

number_of_objects = 100

outputfilename = "1.csv"
filehandle = None

bucketDim = 100000
spatialMap = None

avgN = 10
minN = 3
maxN = 15

allArgs=""
class SpatialMap:
    dim = 0
    buckets = {}
    
    def __init__(self, dim):
        self.dim = dim
    
    def clearBuckets(self):
        for hash in self.buckets:
            self.buckets[hash] = []
    
    # use a tuple as hash
    def getHash(self, loc):
        return (math.floor(loc[0]/self.dim), math.floor(loc[1]/self.dim), math.floor(loc[2]/self.dim))
        
    def addToBuckets(self, co):
        loc = co.getXYZ()
        global CollideScale
        r = co.getScale() * (CollideScale)
        hashes = []
        
        minP = (loc[0]-r, loc[1]-r, loc[2]-r)
        maxP = (loc[0]+r, loc[1]+r, loc[2]+r)
        
        min, max = self.getHash(minP), self.getHash(maxP)
        #printNow(r)
        #printNow(min)
        #printNow(max)
        
        for i in range(min[0], max[0]+1):
            for j in range(min[1], max[1]+1):
                for k in range(min[2], max[2]+1):
                    self.buckets.setdefault( (i, j, k), [] ).append(co)
                    #printNow(len(self.buckets[(i,j,k)]))
                    
        #printNow(self.buckets)
        
    def getBuckets(self, co):
        ret = []
        loc = co.getXYZ()
        global CollideScale
        r = co.getScale() * (CollideScale)
        minP = (loc[0]-r, loc[1]-r, loc[2]-r)
        maxP = (loc[0]+r, loc[1]+r, loc[2]+r)
        
        min, max = self.getHash(minP), self.getHash(maxP)
        
        for i in range(min[0], max[0]+1):
            for j in range(min[1], max[1]+1):
                for k in range(min[2], max[2]+1):
                    if (i,j,k) in self.buckets:
                        ret.append(self.buckets[(i, j, k)])
                    else:
                        sys.stderr.write("Bucket not found!")
                        sys.stderr.flush()
                        assert False
        #printNow(len(ret))
        return ret

def printNow(string, back=False):
    if back:
        sys.stdout.write("\r")
        print(string, end=" ")
    else:
        print(string)
    sys.stdout.flush()

def rrange(min, max):
    return random.random()*(max-min)+min

class CollisionObject:
    mesh = None
    object = None
    vel = (0, 0, 0)
    loc = (0, 0, 0)
    lastloc = (0, 0, 0)
    spawned = False
    hidden = False    
    immune = 0
    mass = 0
    radius = 0
    scale = 0
    envisat = False

    X = (0, 0, 0) #unit vetors for the projected plane
    Y = (0, 0, 0)
    Z = (0, 0, 0)
    P = None
    e = None
    alpha0 = 0
    gamma = 0
    E_M = {None:None}
    a=0

    t = 0
    
    def __init__(self, loc=(None, None, None), mass=None, vel=(None, None, None),envisat=False):
        self.spawn(loc, mass, vel, envisat)
        
    def createOrbitEnvisat(self):
        self.P=6000
        self.e=0
        self.alpha0 = 0
        self.t0 = 0    
        self.a = pow((G*Mearth)*pow(self.P,2)/(4*math.pi*math.pi),(1/3))
        self.createNormal(rrange(-1,1),rrange(-1,1),rrange(-1,1))
        self.envisat = True
        
    def createOrbitR(self):
        self.P=rrange(pLow,pHigh)
        self.e=rrange(0,0.5)
        self.alpha0 = rrange(0,2*math.pi)
        self.t0 = rrange(0,self.P)     
        self.a = pow((G*Mearth)*pow(self.P,2)/(4*math.pi*math.pi),(1/3))
        self.createNormal(rrange(-1,1),rrange(-1,1),rrange(-1,1))
        

        # print("start")
        # print("p",self.P,"e",self.e,"alpha0",self.alpha0,"t0",self.t0,"a",self.a,self.X,self.Y,self.Z)
        
        # self.getCartLoc()
         
        # self.createOrbit()
        # print("end")
        # print("p",self.P,"e",self.e,"alpha0",self.alpha0,"t0",self.t0,"a",self.a,self.X,self.Y,self.Z)

    def createNormal(self,A,B,C):
        self.Z = self.normalize([A,B,C])
        self.X = self.normalize(self.cross([0,1,0],self.Z))
        self.Y = self.normalize(self.cross(self.Z,self.X))
        
    def normalize(self,r):
        len = math.sqrt(pow(r[0],2)+pow(r[1],2)+pow(r[2],2))
        
        return[r[0]/len, r[1]/len, r[2]/len]
    
    
    def cross(self,A,B):
        x = (A[1]*B[2] - A[2]*B[1])
        y = (A[2]*B[0] - A[0]*B[2])
        z = (A[0]*B[1] - A[1]*B[0])
        
        return(x,y,z)
        
    def createOrbit(self):
        v = self.vel
        r = self.loc
        
        A = (r[2]*v[1] - r[1]*v[2])*(r[0]*v[1] - r[1]*v[0])
        B = (r[2]*v[0] - r[0]*v[2])*(r[1]*v[0] - r[0]*v[1])
        C = (r[0]*v[1] - r[1]*v[0])*(r[1]*v[0] - r[0]*v[1])
        

        
        self.createNormal(A,B,C)
        
        #Theta is the clockwise angle from x, Phi is the clockwise angle from y
       
        # self.phi = math.atan2((r[2]*v[1] - r[1]*v[2]),v[0]*r[1]-v[1]*r[0]) + math.pi
        
        # self.phi = math.atan2(r[2] - v[2]*v[1]/r[1],r[0] - v[0]*v[1]/r[1]) 
        
        # self.theta = math.atan2(r[2] - v[2]*v[0]/r[0],r[1] - v[1]*v[0]/r[0]) 
        
        # self.theta = math.atan2(math.cos(self.phi)*r[2] + math.sin(self.phi)*r[0], r[1])

        r2 = pow(r[0],2) + pow(r[1],2) + pow(r[2],2)
        rMag=math.sqrt(r2)
        v2 = pow(v[0],2) + pow(v[1],2) + pow(v[2],2)
        w2 = r2*v2
        
        # print("out",rMag,math.sqrt(v2))

        EperM = v2/2 - G*Mearth/rMag
        
        self.e = math.sqrt(EperM*2*w2/pow(G*Mearth,2) + 1) #eccentricity

        # print(self.e)
        
        if(self.e > 0.9):
            self.despawn()
            return
            #self.e = 0.9
            
        rm = (math.sqrt(pow((G*Mearth),2) + 2*EperM*r2*v2)-G*Mearth)/(2*EperM)

        self.a = rm/(1-self.e)
        
        self.P = math.sqrt(pow(self.a,3)*4*math.pi*math.pi/(G*Mearth)) #Orbital period
    
        #if self.P < pLow or self.P > pHigh:
        #    self.despawn()
        #    return
    
        x = r[0]
        y = r[1]
        z = r[2]
        
        
        arg = ((self.a*(1-pow(self.e,2))/rMag)-1)/self.e
        if(abs(arg) > 1):
            arg = 1
        
        self.alpha0 = math.acos(arg)
           
        E = 2*math.atan(math.sqrt((1-self.e)*pow(math.tan(self.alpha0/2),2)/(1+self.e)))

        self.t0 = (E - self.e*math.sin(E))*2*math.pi/self.P
        
    def Efun(self,E,M):
        return E-self.e*math.sin(E) - M
    
    def getE(self,M):
        M = round(M,8)
        if M in self.E_M :
            return self.E_M[M]
        
        if (self.e > 0.8):
            E = math.pi
        else:
            E = M
 
        while abs(self.Efun(E,M)) > 0.01:
            E = E - (E - self.e*math.sin(E) - M)/(1-self.e*math.cos(E))

        self.E_M[M] = E
        return E
        
    def getOrbitPos(self):
        
        t = (self.t-self.t0)%self.P
        
        M = (t)*2*math.pi/self.P

        E = self.getE(M)
        
        
        alpha = 2*math.atan2(math.sqrt((1+self.e)/(1-self.e))*math.tan(E/2),1)
        dist = self.a*(1-self.e*math.cos(E))
        
        return (dist,alpha)
    
    def getCartLoc(self):

        pos = self.getOrbitPos()
        
        x0 = pos[0]*math.cos(pos[1]-self.alpha0)
        y0 = pos[0]*math.sin(pos[1]-self.alpha0)


        vel2 = G*Mearth*(2/pos[0] - 1/self.a)
        
       

        b = self.a*math.sqrt(1-pow(self.e,2))

        c = math.sqrt(pow(self.a,2) - pow(b,2))
        
        if(y0 ==0):
            vx0 = 0
            vy0 = math.sqrt(vel2)*x0/abs(x0)
        elif(x0 == 0):
            vx0 = -math.sqrt(vel2)*y0/abs(y0)
            vy0 = 0
        else:
            beta = -pow(b/self.a,2)*(x0+c)/y0
            vx0 = -math.sqrt(vel2/(1+pow(beta,2)))*y0/abs(y0)
            vy0 = beta*vx0
   
        # print("in",math.sqrt(pow(vx0,2) + pow(vy0,2)),math.sqrt(pow(x0,2) + pow(y0,2)))
        
        vx = self.X[0]*vx0 + self.Y[0]*vy0
        vy = self.X[1]*vx0 + self.Y[1]*vy0
        vz = self.X[2]*vx0 + self.Y[2]*vy0
        
        self.vel = (vx,vy,vz)
        
     
        x = self.X[0]*x0 + self.Y[0]*y0
        y = self.X[1]*x0 + self.Y[1]*y0
        z = self.X[2]*x0 + self.Y[2]*y0

        self.lastloc = self.loc
        self.loc = (x,y,z)
        
        # print("out",math.sqrt(pow(vx,2) + pow(vy,2) + pow(vz,2)),math.sqrt(pow(x,2) + pow(y,2) + pow(z,2)))
        
    def spawn(self, loc, mass, velocity, envisat):
        self.spawned = True
        self.immune = 0
        
        if envisat:
            self.createOrbitEnvisat()
        elif loc[0] is None or loc[1] is None or loc[2] is None:
            self.createOrbitR()
        else:
            self.loc = loc
            self.vel = velocity
            self.createOrbit()   
            
        if envisat:
            self.setMass(8211)        
        elif mass is None:
            self.setMass(rrange(massMin, massMax))
        else:
            self.setMass(mass)
    
    def hide(self):
        self.hidden = True
    
    def despawn(self):
        self.hide()      
    
    def getR2(self):
        if self.object is not None:
            loc = self.loc
            x = loc[0]
            y = loc[1]
            z = loc[2]
            return (x*x+y*y+z*z)
        return 0
        
    def getXYZ(self):
        loc = self.loc
        return loc
    
    def getlastXYZ(self):
        loc = self.lastloc
        return loc
        
    def getScale(self):
        return self.scale    
             
    def getMass(self):
        return self.mass
        
    def setMass(self, mass):
        self.mass = mass
        self.scale = ((3.*4./math.pi*mass/Density)**(1./3.))
           
    def tick(self):
        if self.hidden:
            return
        if self.immune > 0:
            self.immune -= 1      

        self.t = self.t + DT

        self.getCartLoc()
        global spatialMap
        spatialMap.addToBuckets(self)
        
newcb = []
def checkCollision(a, b, frame):
    if a == b:
        return
    if a.hidden or b.hidden or a.immune > 0 or b.immune > 0:
        return False
        
    global newcb
    
    #printNow("compare")
    
    r1 = a.getXYZ()
    r2 = b.getXYZ()
    
    x1 = r1[0]
    y1 = r1[1]
    z1 = r1[2]
    
    x2 = r2[0]
    y2 = r2[1]
    z2 = r2[2]
    
    #d = min(abs(x2-x1),abs(y2-y1),abs(z2-z1))
    if useinterp:
        ##get the last locations
        r10 = a.getlastXYZ()
        r20 = b.getlastXYZ()

        ##create the line definitions
        #theses are the base point
        l10 = np.array([r10[0], r10[1], r10[2]])
        l20 = np.array([r20[0], r20[1], r20[2]])
        
        #these are the direction vectors
        l1 = np.array([r1[0]-r10[0], r1[1]-r10[1], r1[2]-r10[2]])
        l2 = np.array([r2[0]-r20[0], r2[1]-r20[1], r2[2]-r20[2]])

        #create an arra for utility
        util = np.array([0,0,1])


        #put them all into a matrix as columns, this is used to create the shortest distance vector
        line = np.hstack((l1.T,-l2.T,l10.T-l20.T))

        #put the directions into a matrix as rows, this is used to check perpendicular
        perp = np.vstack((l1,l2,util))

        #this solves for how far along each line the nearest point is
        ab = np.linalg.solve(np.dot(perp, line) , util.T)

        a = ab[0]
        b = ab[1]

        #these should be between 1 and 0 as that is the line segment we are checking

        a = min(max(0,a),1)
        b = min(max(0,b),1)

        #turn these into points
        p = l10 + a*l1
        q = l20 + b*l2

        #get the distance between them
        d = np.numpy.linalg.norm(p-q)        
        
        s1 = a.getScale()
        s2 = b.getScale()
        
        if s1 is None or s2 is None:
            return False
        s = CollideScale*(s1+s2)

        #this is to check if the closest points are within the radius sum
        if d < s: #do the more detailed checking
            #get the lengths of the two lines
            length1 = np.numpy.linalg.norm(l1)
            length2 = np.numpy.linalg.norm(l2)

            #take the average
            lave = (length1 +  length2)/2

            
            #this is to check if the objects were there at the same time. This is approximate, and is most accurate for similar orbits.
            if abs(a-b) < (s/lavg *(1 - d/s)): #if s=d, a=b req'd. if d=0, |a-b| < s/lavg
                collision = True
            else:
                collision = False


        else:
            collision = False


    else:
        d = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        
        s1 = a.getScale()
        s2 = b.getScale()
        
        if s1 is None or s2 is None:
            return False
        s = CollideScale*(s1+s2)
        collision = d < s
    
    if collision:
        #print("boom")
        
        # print(r1,v1)
        # print(r2,v2)
    
        # collision happened.       
        
        # both objects are deleted
        v1 = a.vel
        v2 = b.vel

        vx1 = v1[0]
        vy1 = v1[1]
        vz1 = v1[2]
        
        vx2 = v2[0]
        vy2 = v2[1]
        vz2 = v2[2]       

    
        m1 = a.getMass() 
        m2 = b.getMass()
        
    else: 
        return False
    
    a.hide()
    b.hide()
    
    Mtot = m1 + m2
        
    pxTot = vx1*m1 + vx2*m2
    pyTot = vy1*m1 + vy2*m2
    pzTot = vz1*m1 + vz2*m2
        
        
    # if N is 2, that causes levelling off.
    N = round(random.gauss(avgN, avgN/4))
    if N < minN:
        N = minN
    elif N > maxN:
        N = maxN
    
    #print(N)
    Mtot2 = 0
    
    pxTot2 = 0
    pyTot2 = 0
    pzTot2 = 0
    
    ms = []
    vxs = []
    vys = []
    vzs = []
    
    for i in range(0,N-1):
        ms.append(random.gauss(AverageMass,AverageMass/4))
        if(ms[i] < MinMass):
            ms[i] = MinMass
        vxs.append(random.gauss((vx1 + vx2)/N,(vx1 + vx2)/(4*N)))
        vys.append(random.gauss((vy1 + vy2)/N,(vy1 + vy2)/(4*N)))
        vzs.append(random.gauss((vz1 + vz2)/N,(vz1 + vz2)/(4*N)))
        
    for i in range(0,N-1):
        Mtot2  += ms[i]
        
    for i in range(0,N-1):
        ms[i]  *= Mtot/Mtot2
        
    for i in range(0,N-1):
        pxTot2 += vxs[i]*ms[i]
        pyTot2 += vys[i]*ms[i]
        pzTot2 += vzs[i]*ms[i] 
    
    for i in range(0,N-1):
        vxs[i] *= pxTot/(pxTot2)
        vys[i] *= pyTot/(pyTot2)
        vzs[i] *= pzTot/(pxTot2)

        obj = CollisionObject(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2), ms[i], (vxs[i], vys[i], vzs[i]))
        
        newcb.append(obj)
        
    if Mtot - sum(ms) > 0.01:
        printNow("\n%s -> %s" %(Mtot, sum(ms)))
        assert False

    return True
    
def checkBreakup(a):
    if a.hidden:
        return False
    if not a.envisat or a.t < breakupTime:
        return False
        
    v1 = a.vel

    vx1 = v1[0]
    vy1 = v1[1]
    vz1 = v1[2]
    
    vx2 = 0
    vy2 = 0
    vz2 = 0
    
    m1 = a.getMass() 
    m2 = 0
    
    r1 = a.getXYZ()
        
    x1 = r1[0]
    y1 = r1[1]
    z1 = r1[2]
    
    a.hide()
    
    Mtot = m1
        
    pxTot = vx1*m1 
    pyTot = vy1*m1 
    pzTot = vz1*m1 
        
    N = round(Mtot/AverageMass)
    #print(N)
    Mtot2 = 0
    
    pxTot2 = 0
    pyTot2 = 0
    pzTot2 = 0
    
    ms = []
    vxs = []
    vys = []
    vzs = []
    
    
    for i in range(0,N-1):
        ms.append(random.gauss(AverageMass,AverageMass/4))
        if(ms[i] < MinMass):
            ms[i] = MinMass
        vxs.append(random.gauss((vx1)/N,(vx1)/(4*N)))
        vys.append(random.gauss((vy1)/N,(vy1)/(4*N)))
        vzs.append(random.gauss((vz1)/N,(vz1)/(4*N)))
        
    for i in range(0,N-1):
        Mtot2  += ms[i]
        
    for i in range(0,N-1):
        ms[i]  *= Mtot/Mtot2
        
    for i in range(0,N-1):
        pxTot2 += vxs[i]*ms[i]
        pyTot2 += vys[i]*ms[i]
        pzTot2 += vzs[i]*ms[i] 
    
    for i in range(0,N-1):
        vxs[i] *= pxTot/(pxTot2)
        vys[i] *= pyTot/(pyTot2)
        vzs[i] *= pzTot/(pxTot2)

        obj = CollisionObject(((x1),(y1),(z1)), ms[i], (vxs[i], vys[i], vzs[i]))
        
        newcb.append(obj)
    
    return True

def strTime():
    return time.strftime("%y/%m/%d %H:%m:%S")

def fileWrite(str):
    global outputfilename
    global filehandle
    if filehandle is not None:
        try:
            filehandle.write("%s\n" % str)
        except Exception as e:
            print(str(e))
            
def openFile():
    global outputfilename
    global filehandle
    if filehandle is None:
        try:
            filehandle = open(outputfilename, "a")
        except Exception as e:
            print(str(e))
            filehandle = None
            
def closeFile():
    global filehandle
    if filehandle is not None:
        try:
            filehandle.close()
            filehandle = None
        except Exception as e:
            print(str(e))

last_hidden = 0  
def next_frame():
    global current_frame
    global end_frame
    global newcb
    global cb    
    
    if current_frame == end_frame:
        return False
    
    global spatialMap
    spatialMap.clearBuckets()
    
    for i in range(len(cb)):
        cb[i].tick()
        

    last_p = 0
    newcb = []
    start = time.clock()
    ccount = 0
    bcount = 0
    i = 0
    for A in cb:
        if A.hidden:
            continue
        buckets = spatialMap.getBuckets(A)
        #printNow(len(buckets))
        #printNow(buckets)
        checked = []
        for bucket in buckets:
            for B in bucket:
                if B in checked:
                    continue
                checked.append(B)
                if (checkCollision(A, B, current_frame)):
                    ccount += 1        
        if checkBreakup(A):
            bcount += 1
        p = 100.*i/len(cb)
        if p > last_p + 0.05:
            printNow("[%.2f%%]" % p, back=True)
            last_p = p
        i += 1
             
    printNow("[100.00%]", back=True)
    newcount = len(newcb)
    for i in range(len(newcb)):
        cb.append(newcb[i])
    
    global last_hidden
    hidden = 0
    sizes = []
    for i in range(len(cb)):
        if not cb[i].hidden:
#            sizes.append(str(cb[i].getScale()))
             sizes.append(str(cb[i].getMass()))
        else:
            hidden += 1
    
    new_hidden = hidden - last_hidden - (ccount*2)
    last_hidden = hidden
           
    fileWrite("%f, %d, %d, %d, %s" % (current_frame * DT, len(sizes), ccount, bcount, ' ,'.join(sizes)))
    end = time.clock()
    if newcount > 0:
        pignored = 100.*new_hidden/newcount 
    else:
        pignored = 0.
    printNow("Done step: %d (%.2fs). Num: %d. Collisions: %d. Breakups: %d. New: %d, of which %d were ignored [%.1f%%] (%ds)" % (current_frame, current_frame*DT, len(sizes), ccount, bcount, newcount, new_hidden, pignored, round(end-start)))
    current_frame += 1
    return True

def _start():
    res = True
    while res:
        res = next_frame() 
       
    
def start():
    global outputfilename
    fn = outputfilename
    outputfilename = os.path.join(os.path.dirname(os.path.realpath(__file__)), fn)
    while os.path.exists(outputfilename):
        fn = "_%s" % fn
        outputfilename = os.path.join(os.path.dirname(os.path.realpath(__file__)), fn)
    print(outputfilename)
        
    openFile()
    
    global number_of_objects
    fileWrite("# [%s] Started simulation of %d objects" % (strTime(), number_of_objects))
    fileWrite("# AverageMass: %f" % AverageMass)
    fileWrite("# CollideScale: %f" % CollideScale)
    fileWrite("# Density: %f" % Density)
    fileWrite("# time (s), tot#objects, #collisions, #breakups, masses (kg) ...")
    fileWrite("# %s" % allArgs)
    
    global cb
    for i in range(number_of_objects):
        cb.append(CollisionObject())
    
    global startBreakup
    if(startBreakup):
        obj = CollisionObject(envisat=True)
        cb.append(obj)
    
    global spatialMap
    global bucketDim
    spatialMap = SpatialMap(bucketDim)
    
    _start()
    closeFile()

class HelpOnErrorParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help()
        sys.stderr.write("\n\nERROR: %s\n" % message)
        sys.exit(2)
    
parser = HelpOnErrorParser(description="Kessler Syndrome Simulation.", 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--massMin", "-m", type=float, default=massMin, help="Minimum mass (initial generation) (kg)")
parser.add_argument("--massMax", "-M", type=float, default=massMax, help="Maximum mass (initial generation) (kg)")
parser.add_argument("--collideScale", "-C", type=float, default=CollideScale, help="Collide scale")
parser.add_argument("--averageMass", "-A", type=float, default=AverageMass, help="Average mass (kg)")
parser.add_argument("--minMass", "-i", type=float, default=MinMass, help="Minimum mass (after collision) (kg)")
parser.add_argument("--density", "-d", type=float, default=Density, help="Density (kg/m^3)")
parser.add_argument("--breakup", "-b", type=bool, default=startBreakup, help="Use breakup")
parser.add_argument("--breakupTime", "-B", type=float, default=breakupTime, help="Breakup time (s)")
parser.add_argument("--num", "-N", "-n", type=int, default=number_of_objects, help="Number of Objects", required=True)
parser.add_argument("--out", "-o", default=outputfilename, help="File to output csv data to", required=True)
parser.add_argument("--timestep", "-t", type=float, default=DT, help="Timestep (s)")
parser.add_argument("--end", "--steps", "-e", type=int, default=end_frame, help="Number of time steps to simulate", required=True)
parser.add_argument("--dim", "-D", type=float, default=bucketDim, help="Dimension of spatial segmentation")
parser.add_argument("--avgN", "-a", type=int, default=avgN, help="Average number of objects to create after a collision")
parser.add_argument("--minN", "-j", type=int, default=minN, help="Minimum number of objects to create after a collision")
parser.add_argument("--maxN", "-k", type=int, default=maxN, help="Maximum number of objects to create after a collision")
args = parser.parse_args()    

massMin = args.massMin
massMax = args.massMax
CollideScale = args.collideScale
AverageMass = args.averageMass
MinMass = args.minMass
Density = args.density
startBreakup = args.breakup
breakupTime = args.breakupTime
number_of_objects = args.num
end_frame = args.end
outputfilename = args.out
DT = args.timestep
bucketDim = args.dim
avgN = args.avgN
minN = args.minN
maxN = args.maxN
allArgs = str(args)

start()



