# -*- coding: utf-8 -*-
# -*- encoding=utf-8 -*-

from yade import ymport, utils

#### packing you want to test
PACKING='111_10k'

#### interaction range to adjust for coordination number K
intR=1.245
# 1.245 gives K=10 for 111_10k sample

#### no need to go further down ####

#### material definition
DENS=4000
YOUNG=30e9
ALPHA=0.2
TENS=40e5
COH=40e6
FRICT=18

def sphereMat(): return JCFpmMat(type=1,density=DENS,young=YOUNG,poisson=ALPHA,frictionAngle=radians(FRICT),tensileStrength=TENS,cohesion=COH)

#### import sphere packing
O.bodies.append(ymport.textExt(PACKING+'.spheres',scale=1,material=sphereMat))

R=0
Rmax=0
numSpheres=0.
for o in O.bodies:
 if isinstance(o.shape,Sphere):
  o.shape.color=(0.7,0.5,0.3)
  numSpheres+=1
  R+=o.shape.radius
  if o.shape.radius>Rmax:
   Rmax=o.shape.radius
   
Rmean=R/numSpheres
print 'Rmean=', Rmean

#### engines definition
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom')],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(label='interactionLaw')]
	),
	NewtonIntegrator(damping=0.5),
	
]

###################################################################################
#### default timestep (for interaction creation -> first timestep)
O.dt=0.1*utils.PWaveTimeStep()

#### to manage interaction detection factor during the first timestep and then set default interaction range (intRadius=1)
O.step();

#### coordination number verification and reinforcement of boundary particles
numSSlinks=0
numCohesivelinks=0
numFrictionalLinks=0
for i in O.interactions:
    if not i.isReal : continue
    if isinstance(O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere):
     numSSlinks+=1
     if i.phys.isCohesive :
      numCohesivelinks+=1
     else :
      numFrictionalLinks+=1

print "nbSpheres=", numSpheres," | nblinks=", numSSlinks, " || K =", 2.0*numCohesivelinks/numSpheres

##### YADE interface
#from yade import qt
#v=qt.Controller()
#v=qt.View()

