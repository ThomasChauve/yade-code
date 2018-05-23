# -*- coding: utf-8 -*-
# encoding: utf-8
from yade import ymport, utils, pack, export, plot
from pylab import *
import math

#### SIMULATIONS DEFINED HERE
utils.readParamsFromTable(Y135=16e9,Y24=32e9,SXX=-9.e6,SYY=-6.e6,SZZ=-6e6,FOLDER='test/',NAME='toto',noTableOk=True)
from yade.params.table import *

os.mkdir(FOLDER)


PACK='111_50k'
DFN='penny_R0.1'
intR=1.245

Sxx=SXX # Sigmaxx
Syy=SYY # Sigmayy
Szz=SZZ # Sigmazz

DENS=4000
YOUNG=Y135
YOUNG2=Y24
FRICT=2
RFRICT=0
ALPHA=0.3
TENS=14e6 
COH=32e6


### Fluid properties
KFluid=2.2e11 # bulk modulus of the fluid (1/compressibility)
visc=1 # viscosity of the fluid
pFactor=1.8e-11 # to scale the permeability of the rock matrix: useless if lines 133-136 are not commented (impermeable matrix) -> cf. permeametre.py: 1.8e-11 gives a permeability of 1e-16 m2 for 111_10k
slotAperture=1e-3 # initial aperture of pre-existing fracture where the injection is done

### hydraulic loading
flowRate=1e-5 # injection flow rate

### Simulation Control
saveData=10 # data record interval
iterMax=2e4
saveVTK=10 # number of Vtk files
OUT=FOLDER+PACK+NAME


### preprocessing to get dimensions
O.bodies.append(ymport.text(PACK+'.spheres'))

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
numSpheres=0.
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   numSpheres+=1
   R+=o.shape.radius
   if o.shape.radius>Rmax:
     Rmax=o.shape.radius
Rmean=R/numSpheres

print 'X=',X,' | Y=',Y,' | Z=',Z,' || nbSpheres=',numSpheres,' | Rmean=',Rmean

O.reset()
### material definition
### material definition
def sphereMat(): return JCFpmMat(type=1,density=DENS,young=YOUNG,poisson=ALPHA,tensileStrength=TENS,cohesion=COH,frictionAngle=radians(FRICT),jointNormalStiffness=YOUNG/(pi*Rmean),jointShearStiffness=ALPHA*YOUNG/(pi*Rmean),jointTensileStrength=0.,jointCohesion=0.,jointFrictionAngle=radians(FRICT),jointDilationAngle=radians(0))


def sphereMat2(): return JCFpmMat(type=1,density=DENS,young=YOUNG2,poisson=ALPHA,tensileStrength=TENS,cohesion=COH,frictionAngle=radians(FRICT),jointNormalStiffness=YOUNG2/(pi*Rmean),jointShearStiffness=ALPHA*YOUNG2/(pi*Rmean),jointTensileStrength=0.,jointCohesion=0.,jointFrictionAngle=radians(FRICT),jointDilationAngle=radians(0))



def wallMat(): return JCFpmMat(type=0,density=DENS,young=YOUNG,frictionAngle=radians(0))

### walls
mn,mx=Vector3(xinf+0.1*Rmean,yinf+0.1*Rmean,zinf+0.1*Rmean),Vector3(xsup-0.1*Rmean,ysup-0.1*Rmean,zsup-0.1*Rmean)

walls=utils.aabbWalls(oversizeFactor=1.5,extrema=(mn,mx),thickness=0.1*min(X,Y,Z),material=wallMat)
wallIds=O.bodies.append(walls)

### packing
O.bodies.append(ymport.text(PACK+'.spheres',material=sphereMat))

### DFN
O.bodies.append(ymport.stl(DFN+'.stl',color=(0.9,0.9,0.9),wire=False,material=wallMat))
execfile('identifyInitialFractures.py')

### to reduce initial crack extent: penny crack in the (Y,Z) plane, centre at Yc=0.5 and Zc=0.5, radius 0.1
for o in O.bodies:
 if isinstance(o.shape,Sphere):
   if ( (o.state.pos[2] - 0.5)**2 + (o.state.pos[1] - 0.5)**2 ) > (0.05+Rmean)**2 :
     o.state.onJoint=False

#### engines
### Triaxial Engine
triax=TriaxialStressController(
	internalCompaction=False
	,stressMask=7
	,goal1=Sxx
	,goal2=Syy
	,goal3=Szz
	,max_vel=0.01
)

### Flow Engine
flow=DFNFlowEngine(	
        isActivated=False
        ,useSolver=3
        ,boundaryUseMaxMin = [0,0,0,0,0,0]
        ,bndCondIsPressure = [1,1,1,1,1,1]
        ,bndCondValue=[0,0,0,0,0,0]
        ,permeabilityFactor=pFactor
        ,viscosity=visc
        ,fluidBulkModulus=KFluid
        ### DFN related
        ,clampKValues=False
        ,jointsResidualAperture=slotAperture
        
        ## segfault with following
        #,updatePositions=True
        
        ## segfault with following
        #,defTolerance=-1
        #,updatePositions=True
        #,meshUpdateInterval=1000000000
        
)

### with DFNFlow, we can block every cells not concerned with fractures with the following function: if these lines are commented (permeable matrix), you will get warning about cholmod: is it an issue? I am not sure yet but it is annoying...
def blockStuff():
	for k in range(flow.nCells()): flow.blockCell(k,True)
flow.blockHook="blockStuff()"

### Simulation is defined here (DEM loop, interaction law, servo control, recording, etc...)
O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Box_Aabb(),Bo1_Sphere_Aabb(aabbEnlargeFactor=intR,label='Saabb')]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intR,label='SSgeom'),Ig2_Box_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,label='interactionPhys')],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True,neverErase=1,recordCracks=True,Key=OUT,label='interactionLaw')]
	),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.8,defaultDt=0.1*utils.PWaveTimeStep()),
        triax,
        flow,
        NewtonIntegrator(damping=0.4,label="newton"),
        PyRunner(iterPeriod=int(1),initRun=True,command='crackCheck()',label='check'),
        PyRunner(iterPeriod=int(saveData),initRun=True,command='recorder()',label='recData'),
        PyRunner(iterPeriod=int(1),initRun=True,command='saveFlowVTK()',label='saveFlow',dead=1),
        #PyRunner(iterPeriod=int(1),initRun=True,command='saveAperture()',label='saveAperture',dead=1),
        VTKRecorder(iterPeriod=int(1),initRun=True,fileName=OUT+'-',recorders=['spheres','bstresses','cracks','colors'],Key=OUT,label='saveSolid',dead=0)
]
 
# these lines can be a problem depending on the configuration of your computer
#from yade import qt
#v=qt.Controller()
#v=qt.View()

#### custom functions

#### check if new cracks are created to update "flow mesh permeability"
cks=cks0=0
def crackCheck():
  global tensCks, shearCks, cks, cks0
  cks=interactionLaw.nbTensCracks+interactionLaw.nbShearCracks
  if cks>(cks0): 
    #print 'new crack! Update triangulation!'
    flow.updateTriangulation=True
  cks0=cks

### save flow field (pressure and velocity)
def saveFlowVTK():
 flow.saveVtk(folder='VTK')
 
#### save cracks aperture
#from yade import export
#vtkExporter = export.VTKExporter('cracks')
#def saveAperture():
  #vtkExporter.exportContactPoints(what=[('b','i.phys.isBroken'),('n','i.geom.normal'),('s','i.phys.crossSection'),('a','i.phys.crackJointAperture')])

### save macroscopic data
ex0=ey0=ez0=0
def recorder():
    global ex0,ey0,ez0
    crackVolume=crackSurface=0
    for i in O.interactions:
        if i.phys.isBroken:
	    crackVolume+=i.phys.crossSection*i.phys.crackJointAperture
	    crackSurface+=i.phys.crossSection
    yade.plot.addData( t=O.time
			,i=O.iter
			,ex=triax.strain[0]-ex0
			,ey=triax.strain[1]-ey0
			,ez=triax.strain[2]-ez0
			,sx=triax.stress(triax.wall_right_id)[0]
			,sy=triax.stress(triax.wall_top_id)[1]
			,sz=triax.stress(triax.wall_front_id)[2]
			,p=flow.getPorePressure((xinf+X/2.,yinf+Y/2.,zinf+Z/2.))
			,tc=interactionLaw.nbTensCracks
			,sc=interactionLaw.nbShearCracks
			,p32=crackSurface
			,p33=crackVolume
			,unbF=utils.unbalancedForce() 
    )
    yade.plot.saveDataTxt(OUT)





nbP=0
nb24=0
for o in O.bodies:
	if isinstance(O.bodies[o.id].shape,Sphere):
		nbP+=1
		#print a.bound.refPos[2]
   		if ((O.bodies[o.id].bound.refPos[0]<0.4) and (O.bodies[o.id].bound.refPos[0]>0.2)) or ((O.bodies[o.id].bound.refPos[0]<0.8) and (O.bodies[o.id].bound.refPos[0]>0.6)):
			nb24+=1
			O.bodies[o.id].material=sphereMat2()
			O.bodies[o.id].shape.color=yade.utils.Vector3(1, 0, 0)
		else:
			O.bodies[o.id].shape.color=yade.utils.Vector3(0, 0, 1)

print 'percentage box24=', float(nb24)/nbP

itY135=0
itY24=0
itXX=0
it=0
for o in O.bodies:
	if isinstance(O.bodies[o.id].shape,Sphere):
		it+=1
		if O.bodies[o.id].mat.young==Y135:
			itY135+=1
		elif O.bodies[o.id].mat.young==Y24:
			itY24+=1
		else:
			itXX+=1

print 'percentage Y135 : ', float(itY135)/it ,'| percentage Y24 : ', float(itY24)/it ,'percentage other :', float(itXX)/it
