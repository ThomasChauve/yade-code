# -*- coding: utf-8 -*-
from yade import ymport, utils, pack, export, plot
from pylab import *
import math

#### SIMULATIONS DEFINED HERE

utils.readParamsFromTable(sigma3_conf=-10e6,MAX=1e5, RATE=-0.02,YOUNG=4e9,ALPHA=0.3,TENS=3e6,COH=25e6,OUT='toto',noTableOk=True)

from yade.params.table import *

PACK='111_50k'

PACKING=PACK
OUTPUT=OUT
rate=RATE # loading rate
iterMax=MAX

### Material microproperties
intR=1.245
DENS=4000
finalFricDegree=2


### Mechanical Boundary Conditions
Sxx=sigma3_conf # Sigmaxx
Syy=sigma3_conf # Sigmayy
Szz=sigma3_conf # Sigmazz




### material definition
def sphereMat(): return JCFpmMat(type=1,density=DENS,young=YOUNG,poisson=ALPHA,frictionAngle=radians(finalFricDegree),tensileStrength=TENS,cohesion=COH,jointNormalStiffness=1e9,jointShearStiffness=1e9,jointFrictionAngle=radians(0))

def wallMat(): return JCFpmMat(type=0,density=DENS,young=YOUNG,frictionAngle=radians(0),label='wallMat')

### preprocessing to get dimensions
O.bodies.append(ymport.text(PACKING+'.spheres',scale=1.,shift=Vector3(0,0,0),material=sphereMat))


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
Rmin=0
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

print 'X=',X,' | Y=',Y,' | Z=',Z,' || nbSpheres=',numSpheres,' | Rmean=',Rmean,' | Rmin=',Rmin,' | Rmax=',Rmax




O.reset() 
mn,mx=Vector3(xinf+0.1*Rmean,yinf+0.1*Rmean,zinf+0.1*Rmean),Vector3(xsup-0.1*Rmean,ysup-0.1*Rmean,zsup-0.1*Rmean)
walls=utils.aabbWalls(oversizeFactor=1.5,extrema=(mn,mx),thickness=0.1*min(X,Y,Z),material=wallMat)
wallIds=O.bodies.append(walls)

### packing
O.bodies.append(ymport.text(PACK+'.spheres',material=sphereMat))
print 'density=', sum([b.state.mass for b in O.bodies if isinstance(b.shape,Sphere)])/( (xsup-xinf)* (ysup-yinf)* (zsup-zinf) )
print 'porosity =', ( (xsup-xinf)* (ysup-yinf)* (zsup-zinf) -  4.0/3*pi*sum([b.shape.radius**3 for b in O.bodies if isinstance(b.shape,Sphere)])  ) / ( (xsup-xinf)* (ysup-yinf)* (zsup-zinf) )



triax=TriaxialStressController(
	internalCompaction=True
	,stressMask=7
	,goal1=Sxx
	,goal2=Syy
	,goal3=Szz
	,max_vel=0.01
)


O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom'),Ig2_Box_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True,neverErase=1,recordCracks=True,Key=OUT,label='interactionLaw')]
	),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.8,defaultDt=0.5*utils.PWaveTimeStep()),
        triax,
        NewtonIntegrator(damping=0.25,label="newton")
]

O.dt=0.5*utils.PWaveTimeStep()



  
###first step: compression#######
triax=TriaxialStressController(  
    internalCompaction = True,
    stressMask = 7,
    computeStressStrainInterval = 10,    
    goal1 = Sxx,
    goal2 = Syy,
    goal3 = Szz,
    
)


### Simulation is defined here (DEM loop, interaction law, servo control, recording, etc...)
O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom'),Ig2_Box_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True,neverErase=1,recordCracks=True,Key=OUT,label='interactionLaw')]
	),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.8,defaultDt=0.5*utils.PWaveTimeStep()),
        triax,
        NewtonIntegrator(damping=0.8,label="newton"),
	PyRunner(iterPeriod=int(iterMax/2000),initRun=True,command='recorder()',label='data',dead=True),
	VTKRecorder(iterPeriod=int(iterMax),initRun=True,fileName=OUTPUT+'-',recorders=['spheres','cracks'],Key=OUTPUT,label='vtk')
	
	
]

def recorder():
  e11=-triax.strain[0]  # -
  e22=-triax.strain[1]  # +
  e33=-triax.strain[2]  # -
  s11=-triax.stress(triax.wall_right_id)[0]  # +
  s22=-triax.stress(triax.wall_top_id)[1]  # +
  s33=-triax.stress(triax.wall_front_id)[2]  # +
  unb=unbalancedForce()
  aq=s22-(s11+s33)/2
  aev=-e11-e22-e33
  poro=triax.porosity
  yade.plot.addData({'ae22':e22,
		       'poro':poro,
		       'tc':interactionLaw.nbTensCracks,
		       'sc':interactionLaw.nbShearCracks,
		       'aq':aq,
		       'aev':aev})
  yade.plot.saveDataTxt(OUTPUT)




O.dt=0.5*utils.PWaveTimeStep()


#def confinement_run():
### mechanical loading
while 1:
  O.run(100, True)
  print 'confined state || Sxx=',triax.stress(triax.wall_right_id)[0]/1e6,' | Syy=',triax.stress(triax.wall_top_id)[1]/1e6,' | Szz=',triax.stress(triax.wall_front_id)[2]/1e6,'unbalanced force=',unbalancedForce(), 'poro=',triax.porosity
  if ( unbalancedForce()<0.005*1E-5):# and ((abs(abs(triax.stress(triax.wall_right_id)[0])-abs(Sxx))/abs(Sxx))<0.001) and ((abs(abs(triax.stress(triax.wall_top_id)[1])-abs(Syy))/abs(Syy))<0.001) and ((abs(abs(triax.stress(triax.wall_front_id)[2])-abs(Szz))/abs(Szz))<0.001) ):
    print 'stabilizing || iteration=', O.iter
    O.run(100,True) # to further stabilize the syste
    ex0=triax.strain[0]
    ey0=triax.strain[1]
    ez0=triax.strain[2]   
    break


#def confinement_check_rapid():
print 'confined state || Sxx=',triax.stress(triax.wall_right_id)[0]/1e6,' | Syy=',triax.stress(triax.wall_top_id)[1]/1e6,' | Szz=',triax.stress(triax.wall_front_id)[2]/1e6,'unbalanced force=',unbalancedForce(), 'poro=',triax.porosity
print 'mean velocity=',sum([b.state.vel.norm() for b in O.bodies if isinstance(b.shape,Sphere)])/(len(O.bodies)-6)
print 'mean angvelocity=',sum([b.state.angVel.norm() for b in O.bodies if isinstance(b.shape,Sphere)])/(len(O.bodies)-6)
print 'O.iter=', O.iter




#def deviatoric_load():
##We move to deviatoric loading, let us turn internal compaction off to keep particles sizes constant
triax.internalCompaction=False
## Change contact friction (remember that decreasing it would generate instantaneous instabilities)
#utils.setContactFriction(radians(finalFricDegree))
triax.stressMask = 5
triax.goal2=rate
triax.goal1=Sxx
triax.goal3=Szz


data.dead=False





#confinement_run()
#deviatoric_load()
#deviatoric_write()


O.run(iterMax,True)
