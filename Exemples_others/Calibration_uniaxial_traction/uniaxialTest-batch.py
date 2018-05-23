# -*- coding: utf-8 -*-
from yade import ymport, utils, pack, export, plot
from pylab import *
import math

#### SIMULATIONS DEFINED HERE

utils.readParamsFromTable(PACK='010101_5k',R=1.,Y=27e9,A=0.4,T=6e6,C=15e6,RATE=0.02,OUT='default',MAX=1,noTableOk=True)
from yade.params.table import *

### simulation indices
PACKING=PACK
OUTPUT=OUT

### Simulation Control
rate=RATE # loading rate
### choose how to stop the simulation
iterMax=MAX # maximum number of iteration -> line 144
#maxStress=2.e6 # simulation stops when sigma>maxStress -> lines 146-150

### Material microproperties
intR=R
DENS=4000
YOUNG=Y
ALPHA=A
TENS=T
COH=C
FRICT=0
RFRICT=0

### material definition
def sphereMat(): return JCFpmMat(type=1,density=DENS,young=YOUNG,poisson=ALPHA,frictionAngle=radians(FRICT),tensileStrength=TENS,cohesion=COH,jointNormalStiffness=1e9,jointShearStiffness=1e9,jointFrictionAngle=radians(0))

### preprocessing to get dimensions
O.bodies.append(ymport.text(PACKING+'.spheres',scale=1.,shift=Vector3(0,0,0),material=sphereMat))

dim=utils.aabbExtrema()
xinf=dim[0][0]
X=aabbDim()[0]
yinf=dim[0][1]
Y=aabbDim()[1]
zinf=dim[0][2]
Z=aabbDim()[2]

print 'X=',X,' | Y=',Y,' | Z=',Z

R=0
Rmax=0
numSpheres=0.
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   numSpheres+=1
   R+=o.shape.radius
   if o.shape.radius>Rmax:
     Rmax=o.shape.radius
Rmean=R/numSpheres

### set same color to all spheres -> IMPORTANT FOR REINFORCEMENT TO BE RIGHT!
e=Rmean
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   o.shape.color=(0.7,0.5,0.3)
   # to identify indicators
   if o.state.pos[0]>(xinf+0.5*X-e) and o.state.pos[0]<(xinf+0.5*X+e) and o.state.pos[1]>(yinf+0.75*Y-e) and o.state.pos[1]<(yinf+0.75*Y+e) and o.state.pos[2]>(zinf+0.5*Z-e) and o.state.pos[2]>(zinf+0.5*Z+e):
     yPoint=o.id
     y0=o.state.pos[1]
   if o.state.pos[0]>(xinf+0.75*X-e) and o.state.pos[0]<(xinf+0.75*X+e) and o.state.pos[1]>(yinf+0.5*Y-e) and o.state.pos[1]<(yinf+0.5*Y+e) and o.state.pos[2]>(zinf+0.5*Z-e) and o.state.pos[2]>(zinf+0.5*Z+e):
     xPoint=o.id
     x0=o.state.pos[0]

#print yPoint, y0, xPoint, x0

#### ENGINES DEFINED HERE

### boundary conditions (see utils.uniaxialTestFeatures)
bb=utils.uniaxialTestFeatures(axis=1) # force loading along Y direction
negIds,posIds,axis,crossSectionArea=bb['negIds'],bb['posIds'],bb['axis'],bb['area']

O.engines=[
	ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom')],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(recordCracks=True,Key=OUTPUT,label='interactionLaw')]
	),
	UniaxialStrainer(strainRate=rate,axis=axis,asymmetry=0,posIds=posIds,negIds=negIds,crossSectionArea=crossSectionArea,blockDisplacements=1,blockRotations=1,setSpeeds=0,stopStrain=0.1,label='strainer',dead=True),
	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.5, defaultDt=0.1*utils.PWaveTimeStep()),
	NewtonIntegrator(damping=0.4,label='newton'),
	PyRunner(iterPeriod=int(iterMax/2000),initRun=True,command='recorder()',label='data'),
	VTKRecorder(iterPeriod=int(iterMax),initRun=True,fileName=OUTPUT+'-',recorders=['spheres','cracks'],Key=OUTPUT,label='vtk')
]

epsXX=epsYY=0
def recorder():
    global epsXX, epsYY
    epsXX=(O.bodies[xPoint].state.pos[0]-x0)/(x0-0.5*X)
    epsYY=(O.bodies[yPoint].state.pos[1]-y0)/(y0-0.5*Y)
    yade.plot.addData({'i':O.iter,
		       'eps':strainer.strain,
		       'sigma':strainer.avgStress,
		       'epsXX':epsXX,
		       'epsYY':epsYY,
		       'tc':interactionLaw.nbTensCracks,
		       'sc':interactionLaw.nbShearCracks,
		       'unbF':utils.unbalancedForce()})
    yade.plot.saveDataTxt(OUTPUT)

#### SIMULATION STARTS HERE

### manage interaction detection factor during the first timestep and then set default interaction range
O.step();
### initializes the interaction detection factor
SSgeom.interactionDetectionFactor=-1.
Saabb.aabbEnlargeFactor=-1.

### reinforce bonds for boundary particles coordination number verification
layerSize=0.2
for o in O.bodies:
  if isinstance(o.shape,Sphere):
    if ( o.state.pos[axis]<(dim[0][axis]+layerSize*(dim[1][axis]-dim[0][axis])) ) or ( o.state.pos[axis]>(dim[1][axis]-layerSize*(dim[1][axis]-dim[0][axis])) ) :
      o.shape.color=(1,1,1)
      
### coordination number verification
numSSlinks=0
numCohesivelinks=0
for i in O.interactions:
    if isinstance(O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere):
      numSSlinks+=1
      if ( rate>0 and (O.bodies[i.id1].shape.color==(1,1,1) or O.bodies[i.id2].shape.color==(1,1,1)) ) :
	i.phys.FnMax*=100
	i.phys.FsMax*=100
    if i.phys.isCohesive :
      numCohesivelinks+=1

print "nbSpheres=", numSpheres," | Rmean=", Rmean, " || Kcohesive=", 2.0*numCohesivelinks/numSpheres

### simulation really starts here
strainer.dead=False

O.run(iterMax,True)
