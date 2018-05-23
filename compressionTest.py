# -*- coding: utf-8 -*-
from yade import ymport, utils, pack, export, plot
from pylab import *
import math

#### SIMULATIONS DEFINED HERE

PACK='12.51_20k'
beddingAngle=45 # with respect to the XY plane (perpendicular to the loading here)
R=1.234 # to ensure K=10 (avrg nb of bonds per particle) for this particular sample
RATE=-0.1 # loading rate
MAXI=1e5 # iteration max
VTK=5 # nb of vtk files to record
OUT='12.51_10k_g0_compression_r0.1'

#### Simulation Control
rate=RATE
maxIter=MAXI
vtkSave=int(maxIter/VTK)
PACKING=PACK
intR=R
OUTPUT=OUT

#### Material microproperties (Tournemire argilite)
DENS=4000
YOUNG=16e9
FRICT=2
RFRICT=0
ALPHA=0.3
TENS=14e6 
COH=32e6

## weak planes for TI materials -> CAREFUL: SET smoothJoint=TRUE in Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM
GAMMA=beddingAngle  # angle of bedding with respect to the horizontal direction (0 if planes are perpendicular to loading (horizontal))
DGAMMA=55 # tolerance around GAMMA / cone interval for selecting contacts
BEDNSTIFF=0.2 # normal stiffness of bedding
BEDSSTIFF=1*BEDNSTIFF # shear stiffness of bedding
BEDTENS=0 # tensile strength on bedding (put 1 instead of 0 if you wan to test 0 -> it will let you see the contacts with Paraview)
BEDCOH=0 # cohesion on bedding (put 1 instead of 0 if you wan to test 0 -> it will let you see the contacts with Paraview)
BEDFRIC=0 # friction angle on bedding
BEDDIL=0 # dilation angle on bedding

#### material definition
def sphereMat(): return JCFpmMat(type=1,density=DENS,young=YOUNG,poisson=ALPHA,frictionAngle=radians(FRICT),residualFrictionAngle=radians(RFRICT),tensileStrength=TENS,cohesion=COH)

#### import packing
O.bodies.append(ymport.text(PACKING+'.spheres',scale=1.,shift=Vector3(0,0,0),material=sphereMat))

dim=utils.aabbExtrema()
xmin=dim[0][0]
xmax=dim[1][0]
X=aabbDim()[0]
ymin=dim[0][1]
ymax=dim[1][1]
Y=aabbDim()[1]
zmin=dim[0][2]
zmax=dim[1][2]
Z=aabbDim()[2]

#### set same color to all spheres
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   o.shape.color=(0.7,0.5,0.3)

#### ENGINES DEFINED HERE

## boundary conditions (see utils.uniaxialTestFeatures)
bb=utils.uniaxialTestFeatures()
negIds,posIds,axis,crossSectionArea=bb['negIds'],bb['posIds'],bb['axis'],bb['area']

O.engines=[
	ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb'),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom')],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(recordCracks=True,smoothJoint=True,Key=OUTPUT,label='interactionLaw')]
	),
	UniaxialStrainer(strainRate=rate,axis=axis,asymmetry=0,posIds=posIds,negIds=negIds,crossSectionArea=crossSectionArea,blockDisplacements=1,blockRotations=1,setSpeeds=0,stopStrain=0.1,label='strainer',dead=True),
	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.5, defaultDt=0.1*utils.PWaveTimeStep()),
	NewtonIntegrator(damping=0.4,label='newton'),
	PyRunner(iterPeriod=100,initRun=True,command='recorder()',label='data'),
	PyRunner(iterPeriod=int(1),initRun=False,command='position()',label='pos'),
        VTKRecorder(iterPeriod=int(1),initRun=False,fileName=OUTPUT+'-',recorders=['spheres','colors','bstresses','cracks'],Key=OUTPUT,label='vtk')
]

def recorder():
    yade.plot.addData({'i':O.iter,
		       'eps':strainer.strain,
		       'sigma':strainer.avgStress,
		       'tc':interactionLaw.nbTensCracks,
		       'sc':interactionLaw.nbShearCracks,
		       'tce':interactionLaw.totalTensCracksE,
		       'sce':interactionLaw.totalShearCracksE,
		       'cSurf':interactionLaw.totalCracksSurface,
		       'unbF':utils.unbalancedForce()})
    yade.plot.saveDataTxt(OUTPUT)
    
def position():
  export.text(OUT+'_'+str(O.iter))

#### Initialization of contacts (setting the bonds between the particles and modifying them accordingly to the anisotropy)

### manage interaction detection factor during the first timestep and then set default interaction range
O.step();
### initializes the interaction detection factor
SSgeom.interactionDetectionFactor=-1.
Saabb.aabbEnlargeFactor=-1.

### Introduce microfractures in the packing -> CAREFUL: SET smoothJoint=TRUE
BEDNORM=Vector3(sin(radians(GAMMA)),cos(radians(GAMMA)),0) # normal to bedding (vertical if planes are horizontal)
print 'angle of bedding with respect to horizontal =', GAMMA, ' | normal =', BEDNORM, ' | angle range=', DGAMMA
## search for preferentially oriented contacts and assignation of bedding properties
nbInteractions=0
nbMicroFrac=0
cumLengthMicroFrac=0
cumSurfaceMicroFrac=0
for i in O.interactions:
  nbInteractions+=1
  if i.phys.isOnJoint==False and (O.bodies[i.id1].shape,Sphere) and isinstance(O.bodies[i.id2].shape,Sphere):
    pdct=(i.geom.normal).dot(BEDNORM)
    if abs(pdct) > cos(radians(DGAMMA)): # condition that ensures the contact normal belongs to the cone defined by bedding's normal with more or less dGamma degrees
      ## contact reorientation -> use of smooth contact model
      # bodies
      O.bodies[i.id1].state.onJoint=True
      O.bodies[i.id2].state.onJoint=True
      O.bodies[i.id1].state.joint=1
      O.bodies[i.id2].state.joint=1
      O.bodies[i.id1].state.jointNormal1=BEDNORM
      O.bodies[i.id2].state.jointNormal1=BEDNORM      
      O.bodies[i.id1].mat.jointNormalStiffness=BEDNSTIFF*i.phys.kn/i.phys.crossSection
      O.bodies[i.id2].mat.jointNormalStiffness=BEDNSTIFF*i.phys.kn/i.phys.crossSection
      O.bodies[i.id1].mat.jointShearStiffness=BEDSSTIFF*i.phys.ks/i.phys.crossSection
      O.bodies[i.id2].mat.jointShearStiffness=BEDSSTIFF*i.phys.ks/i.phys.crossSection
      O.bodies[i.id1].mat.jointFrictionAngle=radians(BEDFRIC)
      O.bodies[i.id2].mat.jointFrictionAngle=radians(BEDFRIC)
      O.bodies[i.id1].mat.jointDilationAngle=radians(BEDDIL)
      O.bodies[i.id2].mat.jointDilationAngle=radians(BEDDIL)
      O.bodies[i.id1].mat.jointTensileStrength=BEDTENS*i.phys.FnMax/i.phys.crossSection
      O.bodies[i.id2].mat.jointTensileStrength=BEDTENS*i.phys.FnMax/i.phys.crossSection
      O.bodies[i.id1].mat.jointCohesion=BEDCOH*i.phys.FsMax/i.phys.crossSection
      O.bodies[i.id2].mat.jointCohesion=BEDCOH*i.phys.FsMax/i.phys.crossSection
      # interactions
      i.phys.isOnJoint=True
      i.phys.jointNormal=BEDNORM
      i.phys.jointNormal=BEDNORM*sign(i.phys.jointNormal.dot(i.geom.normal))
      i.phys.initD = abs((O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos).dot(i.phys.jointNormal))
      i.phys.kn*=BEDNSTIFF
      i.phys.ks*=BEDSSTIFF
      i.phys.tanFrictionAngle=tan(radians(BEDFRIC))
      i.phys.tanDilationAngle=tan(radians(BEDDIL))
      i.phys.FnMax*=BEDTENS
      i.phys.FsMax*=BEDCOH
      # statistics
      nbMicroFrac+=1
      cumLengthMicroFrac+=2*sqrt(i.phys.crossSection/pi)
      cumSurfaceMicroFrac+=i.phys.crossSection
      ## comment the following lines if you want to see cracks form anisotropy search
      #if BEDCOH==0 or BEDTENS==0:
	#i.phys.isCohesive=0

print 'nbInteractions=', nbInteractions,' | nbMicroFrac=',nbMicroFrac, ' | cumLengthMicroFrac=', cumLengthMicroFrac, '| cumSurfaceMicroFrac=',cumSurfaceMicroFrac

#### simulation starts here

# loading
O.step();
strainer.dead=False
vtk.iterPeriod=vtkSave
pos.iterPeriod=vtkSave

O.run(int(maxIter),True)