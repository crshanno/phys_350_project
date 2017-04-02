import bpy
import random
import math
import time
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy.optimize import fsolve

G = 6.67408*math.pow(10, -20) #km^3kg^-1s^-2 (Gravitational constant)
Mearth = 5.972*math.pow(10, 24) #kg (mass of earth)

Re = 6371 #km (radius of earth)
framerate = 24 

random.seed(time.time())

vstart=2
rstart = 1

pLow = 88*60
pHigh = 127*60

def rrange(min, max):
    return random.random()*(max-min)+min

class CollisionObject:
    mesh = None
    object = None
    velocity = (0, 0, 0)
    loc = (0, 0, 0)
    spawned = False
    hidden = False    
    immune = 0

    theta = 0
    phi = 0
    P = None
    e = None
    alpha0 = 0
    gamma = 0
    E_M = {None:None}
    a=0

    t = 0
    
    def __init__(self, loc=(None, None, None), scale=None, velocity=(None, None, None)):
        self.spawn(loc, scale, velocity)

    def createOrbitR(self):
        self.P=rrange(pLow,pHigh)
        self.e=rrange(0,0.3)
        self.alpha0 = rrange(0,2*math.pi)
        self.theta = rrange(0,math.pi)
        self.phi = rrange(0,2*math.pi)
        self.t0 = rrange(0,self.P)
        
        self.a = pow((G*Mearth)*pow(self.P,2)/(4*math.pi*math.pi),(1/3))

    def createOrbit(self):
        v = self.velocity
        r = self.loc
        #These are the angles of the normal vector in polar coordinates
        self.phi = math.atan((v[0]*v[2] - r[0]*r[2])/(r[1]*r[2] - v[1]*v[2])) 
        self.theta = math.atan(math.sqrt(pow((v[1]*v[2] -r[1]*r[2]),2) + pow((v[0]*v[2] -r[0]*r[2]),2))/(r[0]*r[1] - v[0]*v[1]))

        r2 = r[0]^2 + r[1]^2 + r[2]^2
        rMag=math.sqrt(r2)
        v2 = pow(v[0],2) + pow(v[1],2) + pow(v[2],2)
        w2 = r2*v2
    
        self.e = math.sqrt(w2*(v2/(G*Mearth) - 2/rMag) + 1) #eccentricity

        if(self.e > 1):
            self.despawn()
            return
        
        EperM = v2/2 - G*Mearth/rMag
        
        rm = (math.sqrt(pow((G*Mearth),2) + 2*EperM*r2*v2)-G*Mearth)/(2*EperM)

        self.a = rm/(1-self.e)
        
        self.P = math.sqrt(pow(self.a,3)*4*math.pi*math.pi/(G*Mearth)) #Orbital period

        self.alpha0 = math.acos((rMag*v2/(G*Mearth)-1)/self.e)-math.atan(r[1]/r[0]) #Angle of r_min from x axis

        theta = math.acos((rMag*v2/(G*Mearth)-1)/self.e)
        
        E = 2*math.atan(math.sqrt((1-self.e)*pow(math.tan(theta/2),2)/(1+self.e)))

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
        
        print(E)
        return (dist,alpha)
    
    def getCartLoc(self):

        pos = self.getOrbitPos()
        
#        vel = sqrt(G*Mearth*(2/pos[0] - 1/self.a))
        
#        rpert = r*self.a*self.e*math.sin(alpha)(1-self.e^2)/pow((1 + self.e*math.cos(alpha)),2)

        x0 = pos[0]*math.cos(pos[1]-self.alpha0)
        y0 = pos[0]*math.sin(pos[1]-self.alpha0)

        x = x0*math.cos(self.theta)*math.cos(self.phi) - y0*math.cos(self.theta)*math.sin(self.phi)
        y = x0*math.cos(self.theta)*math.sin(self.phi) + y0*math.cos(self.theta)*math.cos(self.phi)
        z = -x0*math.sin(self.theta)*math.cos(self.phi) -y0*math.sin(self.theta)*math.sin(self.phi)
        return (x,y,z)
  
    def pickStartingVelocity(self):
        x = rrange(-vstart, vstart)
        y = rrange(-vstart, vstart)
        z = rrange(-vstart, vstart)
        self.velocity = (x, y, z)
    
    def pickStartingLocation(self):
        x = rrange(-rstart, rstart)
        y = rrange(-rstart, rstart)
        z = rrange(-rstart, rstart)
        return (x, y, z)
        
    def spawn(self, loc, scale, velocity):
        self.spawned = True
        self.immune = 10
        self.mesh = bpy.data.meshes.new_from_object(
            bpy.context.scene, 
            bpy.data.objects["CD"],
            True,
            'PREVIEW')
        self.object = bpy.data.objects.new("Object", self.mesh)
        bpy.context.scene.objects.link(self.object)
        if loc[0] is None or loc[1] is None or loc[2] is None:
            self.setLocation(self.pickStartingLocation())
        else:
            x = loc[0]
            y = loc[1]
            z = loc[2]
            self.setLocation((x, y, z))       
        if scale is None:
            self.setScale(rrange(ScaleMin, ScaleMax))
        else:
            self.setScale(scale)
        if velocity[0] is None or velocity[1] is None or velocity[2] is None:
            self.pickStartingVelocity()
        else:
            self.velocity = velocity

        if self.P is None or self.e is None:
            self.createOrbitR()

        self.object.hide = True
        self.object.hide_render = True
        self.object.keyframe_insert(data_path="hide", frame=bpy.context.scene.frame_current-1)  
        self.object.keyframe_insert(data_path="hide_render", frame=bpy.context.scene.frame_current-1)         
        self.object.hide = False
        self.object.hide_render = False
        self.object.keyframe_insert(data_path="hide")
        self.object.keyframe_insert(data_path="hide_render")
    
    def hide(self):
        self.hidden = True
        self.object.hide = True
        self.object.hide_render = True
        self.object.keyframe_insert(data_path="hide")
        self.object.keyframe_insert(data_path="hide_render")  
    
    def despawn(self):
        if self.spawned and self.object is not None:
            self.spawned = False
            self.object.hide = False
            self.object.select = True
            bpy.ops.object.delete()
    
    def getR2(self):
        if self.object is not None:
            loc = self.loc
            x = loc[0]
            y = loc[1]
            z = loc[2]
            return (x*x+y*y+z*z)
        return 0
        
    def getXYZ(self):
        loc = self.object.location
        return loc
        
    def setLocation(self, loc):
        if self.object is not None:
            loc = (loc[0]/100, loc[1]/100, loc[2]/100)
            self.object.location = loc
            self.object.keyframe_insert(data_path="location")            
    
    def getScale(self):
        if self.object is not None:
            return self.object.scale[1]
            
    def setScale(self, scale):
        if self.object is not None:
            self.object.scale = (scale, scale, scale)
    
    def getMass(self):
        scale = self.object.scale[1]
        return density*scale*scale*scale
    
    def _tick(self, rewind=False):
        if self.hidden:
            return
        if rewind:
            print("Rewind not implemented")
            return
        if self.immune > 0:
            self.immune -= 1
            
        dt = 100

        self.t = self.t+dt
        
        self.setLocation(self.getCartLoc())
        
                
    def tick(self, ticks):
        if not self.spawned:
            return
        if ticks == 0:
            pass
        elif ticks < 0:
            print("Rewind not implemented")
        else:
            for i in range(ticks):
                self._tick()

cb = []    
allcb = []
last_frame = 0

def frame_change_handler(scene):
    global last_frame
    frame = scene.frame_current
    if frame == scene.frame_end:
        bpy.ops.screen.animation_cancel(restore_frame=False)
        bpy.app.handlers.frame_change_pre.clear()
        global started
        started = False
        return
    dframe = frame-last_frame
    
    for i in range(len(cb)):
        cb[i].tick(dframe)
        
    deli = []
    delk = []
    # very inefficient N^N collision checking
    for i in range(len(cb)):
       for k in range(len(cb)-i):
          if i != k:
             checkCollision(cb[i], cb[k])
                        
    last_frame = frame

def checkCollision(a, b):
    if a.hidden or b.hidden or a.immune > 0 or b.immune > 0:
        return False
    
    r1 = a.getXYZ()
    r2 = b.getXYZ()
    
    x1 = r1[0]
    y1 = r1[1]
    z1 = r1[2]
    
    x2 = r2[0]
    y2 = r2[1]
    z2 = r2[2]
    
    vr1 = a.velocity[0]
    vtheta1 = a.velocity[1]
    vphi1 = a.velocity[2]
        
    vr2 = b.velocity[0]
    vtheta2 = b.velocity[1]
    vphi2 = b.velocity[2]
    
    d = math.pow(x2-x1,2) + math.pow(y2-y1,2) + math.pow(z2-z1,2)
    d = math.sqrt(d)
    
    s1 = a.getScale()
    s2 = b.getScale()
    if s1 is None or s2 is None:
        return False
    
    s = (s1+s2)/2
    
    if d < s:
        # collision happened.       
        
        # both objects are deleted
        a.hide()
        b.hide()
        
        # create 4 objects. very fake
        a = CollisionObject(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2), (s1+s2)/4, (-vr1, vtheta1, vphi1))
        b = CollisionObject(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2), (s1+s2)/4, (vr1, vtheta1, vphi1))
        c = CollisionObject(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2), (s1+s2)/4, (vr2, -vtheta2, vphi2))
        d = CollisionObject(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2), (s1+s2)/4, (vr2, vtheta2, -vphi2))
        cb.append(a)
        cb.append(b)
        cb.append(c)
        cb.append(d)
        allcb.append(a)
        allcb.append(b)
        allcb.append(c)
        allcb.append(d)
            
        return True
    return False

cb = []    
last_frame = 0

started = False

      
Rmin = 10
Rmax = 20
ScaleMin = 0.1
ScaleMax = 1
Me=1
density = 100  

class KesslerSyndromeStart(bpy.types.Operator):
    bl_idname = "ks.start"
    bl_label = "Start"
    
    def invoke(self, context, event):
        global last_frame
        global cb
        global allcb
        global started
        
        if started:
            return {"FINISHED"}
        started = True
        
        global Rmin
        Rmin = context.scene.Rmin
        global Rmax
        Rmax = context.scene.Rmax
        global ScaleMin
        ScaleMin = context.scene.ScaleMin
        global ScaleMax
        ScaleMax = context.scene.ScaleMax
        global Me
        Me = context.scene.Me
        global density
        density = context.scene.density        
        
        bpy.ops.screen.frame_jump(1)
        
        for i in range(len(allcb)):
            allcb[i].despawn()        
        allcb = []
        cb = []    
        
        for i in range(context.scene.number_of_objects):
            cb.append(CollisionObject())
            allcb.append(cb[i])
            
        bpy.app.handlers.frame_change_pre.clear()
        bpy.app.handlers.frame_change_pre.append(frame_change_handler)

        
        bpy.ops.screen.animation_play()
        return {"FINISHED"}
    
class KesslerSyndromeContinue(bpy.types.Operator):
    bl_idname = "ks.continue"
    bl_label = "Continue"
    
    def invoke(self, context, event):
        bpy.app.handlers.frame_change_pre.clear()
        bpy.app.handlers.frame_change_pre.append(frame_change_handler)
        bpy.ops.screen.animation_play()
        return {"FINISHED"}
        
class KesslerSyndromeStop(bpy.types.Operator):
    bl_idname = "ks.stop"
    bl_label = "Stop"
    
    def invoke(self, context, event):
        global started
        started = False
        bpy.app.handlers.frame_change_pre.clear()
        bpy.ops.screen.animation_cancel(restore_frame=False)
        return {"FINISHED"}
    
class KesslerSyndromeClear(bpy.types.Operator):
    bl_idname = "ks.clear"
    bl_label = "Clear"
    
    def invoke(self, context, event):
        global cb
        global allcb
        for i in range(len(allcb)):
            allcb[i].despawn()
        cb = []
        allcb = []
        return {"FINISHED"}

class KSPanel(bpy.types.Panel):
    bl_label = "Kessler Syndrome"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"
 
    def draw(self, context):
        col1 = self.layout.column(align=True)
        col1.prop(context.scene, "number_of_objects", slider=True)
        r1 = col1.row(align=True)
        r1.prop(context.scene, "Rmin", slider=True)
        r1.prop(context.scene, "Rmax", slider=True)
        r2 = col1.row(align=True)       
        r2.prop(context.scene, "ScaleMin", slider=True)
        r2.prop(context.scene, "ScaleMax", slider=True)
        col1.prop(context.scene, "Me", slider=True)
        col1.prop(context.scene, "density", slider=True)
        col2 = self.layout.column(align=True)
        col2.operator("ks.start")
        col2.operator("ks.continue")
        col2.operator("ks.stop")
        col2.operator("ks.clear")   

def register():
    bpy.utils.register_class(KSPanel)
    bpy.utils.register_class(KesslerSyndromeStart)
    bpy.utils.register_class(KesslerSyndromeContinue)
    bpy.utils.register_class(KesslerSyndromeStop)
    bpy.utils.register_class(KesslerSyndromeClear)
    bpy.types.Scene.number_of_objects = bpy.props.IntProperty \
      (
        name = "N",
        description = "Number of objects to generate",
        default = 10,
        min = 0
      )      
    bpy.types.Scene.Rmin = bpy.props.FloatProperty \
      (
        name = "Rmin",
        description = "Minimum radius to generate objects",
        default = 10,
        min = 0
      )
    bpy.types.Scene.Rmax = bpy.props.FloatProperty \
      (
        name = "Rmax",
        description = "Maximum radius to generate objects",
        default = 20,
        min = 0
      )
    bpy.types.Scene.ScaleMin = bpy.props.FloatProperty \
      (
        name = "ScaleMin",
        description = "Minimum scale to generate objects",
        default = 0.1,
        min = 0
      )
    bpy.types.Scene.ScaleMax = bpy.props.FloatProperty \
      (
        name = "ScaleMax",
        description = "Maximum scale to generate objects",
        default = 1,
        min = 0
      )
    bpy.types.Scene.Me = bpy.props.FloatProperty \
      (
        name = "Me",
        description = "Mass of earth",
        default = 1, #5.972*math.pow(10, 24)
        min = 0
      )
    bpy.types.Scene.density = bpy.props.FloatProperty \
      (
        name = "density",
        description = "Density of objects",
        default = 100,
        min = 0
      )
      
def unregister():
    bpy.utils.unregister_class(KSPanel)
    bpy.utils.unregister_class(KesslerSyndromeStart)
    bpy.utils.unregister_class(KesslerSyndromeContinue)
    bpy.utils.unregister_class(KesslerSyndromeStop)
    bpy.utils.unregister_class(KesslerSyndromeClear)
    
if __name__ == "__main__" :
    register()
