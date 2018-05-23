# -*- coding: utf-8 -*-

from yade import pack, ymport

#### packing you want to test
PACK='111_10k'
intR=1.245 # see coordinationNumber.py

#### matrix permeability: allow you to scale the permeability of the medium (if pFactor=1, no size scaling)
pFactor=1.80e-11 
# 1.80e-11 gives a macroscopic permeability k=1e-16 [m2] for 111_10k

#### material definition
def sphereMat(): return JCFpmMat(type=1,density=4000,young=1e9,poisson=0.3,frictionAngle=radians(30),tensileStrength=1e6,cohesion=1e6)
def wallMat(): return JCFpmMat(type=0,density=4000,young=1e9,poisson=0.3,frictionAngle=radians(0))

#### preprocessing to get dimensions
O.bodies.append(ymport.text(PACK+'.spheres',material=sphereMat))

dim=utils.aabbExtrema()
xinf=dim[0][0]
xsup=dim[1][0]
X=xsup-xinf
yinf=dim[0][1]
ysup=dim[1][1]
Y=ysup-yinf
zinf=dim[0][2]
zsup=dim[1][2]
Z=zsup-zinf

R=0
Rmax=0
Rmin=1000
numSpheres=0.
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   numSpheres+=1
   R+=o.shape.radius
   if o.shape.radius>Rmax:
     Rmax=o.shape.radius
   if o.shape.radius<Rmin:
     Rmin=o.shape.radius
Rmean=R/numSpheres

print "nb spheres=",numSpheres, " | mean Diameter=", 2*Rmean, ",X,Y,Z,", X,Y,Z

O.reset()

### simulation is built up now
mn,mx=Vector3(xinf,yinf,zinf),Vector3(xsup,ysup,zsup)
walls=aabbWalls([mn,mx],oversizeFactor=1.5,thickness=min(X,Y,Z)/100.,material=wallMat)
wallIds=O.bodies.append(walls)

O.bodies.append(ymport.text(PACK+'.spheres',material=sphereMat))

flow=FlowEngine(
        isActivated=1
        ,useSolver=3 # 3 should be used by default as it is the most efficient
        ,boundaryUseMaxMin = [0,0,0,0,0,0]
        ,permeabilityFactor=pFactor
        ,viscosity=0.001
)

O.engines=[
	ForceResetter()
	,InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR)])
	,InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR),Ig2_Box_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1)],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True)]
	)
	,flow
	,GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.8,defaultDt=0.1*utils.PWaveTimeStep())
	,NewtonIntegrator(damping=1)
	,VTKRecorder(iterPeriod=1,initRun=True,fileName='permeametre-',recorders=['spheres','intr'])
]

## if you want to block the particles (better for permeability test)
for o in O.bodies:
 o.state.blockedDOFs+='xyzXYZ'
 
## flow along X direction
print 'Flow along X!'
flow.bndCondIsPressure = [1,1,0,0,0,0]
flow.bndCondValue = [1,0,0,0,0,0]
O.step()
Qin = flow.getBoundaryFlux(0) 
Qout = flow.getBoundaryFlux(1)
#### getFlux gives total discharge -> Q in (m3/s) and we can compute permeability k in [m2] as k=Q*mu*Length/(Area*deltaP)
permeability = abs(Qout)*flow.viscosity*X/(Y*Z) # !!! works if Pout=1 and Pin=0 (deltaP=1)
# Rk: we can also compute k as: k=flow.averageVelocity()*flow.viscosity*X # !!! works if Pout=1 and Pin=0 (deltaP=1)
# conductivity K in [m/s] can be computed as K=rho*g*k/nu
conductivity = permeability*1000*9.82/flow.viscosity 
print "Qin=",Qin," Qout=",Qout," ARE THEY EQUAL? "
print "Permeability [m2]=",permeability," || Hydraulic conductivity [m/s]=",conductivity
flow.saveVtk() # if you want to see the result in Paraview
