#=================================================================================================
# In silico actuation performance investigation of dielectric elastomers with TPMS geometries
# Link: https://doi.org/10.1016/j.euromechsol.2024.105540
#=================================================================================================

# Part of the code for creating TPMS structures was developed by Fayyaz Nosouhi (dehnavifn@gmail.com) and Saeed Khaleghi (saeedkhaleghi123@gmail.com).
# The remaining parts of the code, including the definition of DEA composites and their implementation in ABAQUS, were developed by Mohammad Ali Safaei (mohammadsf1998@gmail.com).
# This Python script requires the use of a UEL subroutine, which implements the constitutive equations for a DEA element. In this regard, the UEL subroutine developed by Ehsan Hajiesmaili was utilized. 
# University of Tehran, December 2024

#=================================================================================

from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import numpy as np
import time
import math
from operator import itemgetter 
import re
import numpy as np

#==========================================================
# 				     File Execution
#==========================================================

execfile('vol.py') # vol.py contains lists that show the corresponding volume fraction for each t (surface level parameter).
def MyRVE(Nxyz,cxyz, RVE_Type, Vol):
	#--------------------------------------------------
	# This function creates a 3D binary voxel structure(0 for void and 1 for solid).
	# The inputs are:
	# Nxyz: number of voxels in x-direction
	
	# cx: number of cells in x-direction
	
	# RVE_Type: the TPMS structure. it can be one of these: 'IWP', 'Schwartz', 'Gyroid', 'Diamond', 'Neovius', 'FKS' or 'FRD'
	# Vol: Target Volume Fraction of the Solid Phase.
	#--------------------------------------------------
	xyz=cxyz*cxyz*cxyz*Nxyz*Nxyz*Nxyz # Total number of voxels
	Lx=2*cxyz*pi # Length of the RVE (x-direction)
	Ly=2*cxyz*pi # Length of the RVE (y-direction)
	Lz=2*cxyz*pi # Length of the RVE (z-direction)
	c1=(Lx)/(cxyz*Nxyz) # Length of each Voxel in x-direction
	c2=(Ly)/(Nxyz*cxyz) # Length of each Voxel in y-direction
	c3=(Lz)/(cxyz*Nxyz) # Length of each Voxel in z-direction
	# Select the appropriate [VolumeFraction-t] list based on the type of structure
	if   (RVE_Type== 'IWP'):
		vol_t = iwp_vol
	elif (RVE_Type== 'Schwarz'):
		vol_t = sch_vol
	elif (RVE_Type== 'Octo'):
		vol_t = oc_vol
	elif (RVE_Type== 'Gyroid'):
		vol_t = gyr_vol
	elif (RVE_Type== 'Diamond'): 
		vol_t = dim_vol
	elif (RVE_Type== 'F'):
		vol_t = f_vol
	elif (RVE_Type== 'Neovius'):  
		vol_t = neov_vol
	elif (RVE_Type== 'CD'):
		vol_t = cd_vol
	elif (RVE_Type== 'FRD'):  
		vol_t = frd_vol
	else: #FKS
		vol_t = fks_vol
	# Find t corresponding to the volume fraction	
	t=0
	Q=0
	for q in range (int(len(vol_t)/2)):
		if vol_t[q*2+1]>Vol:
			t=vol_t[2*q-2]
			Q=q-1
			break
	# Create the Voxel Model
	# if coordinate of the center of a voxel satisfy "F(x,y,z)<t", that voxel will be solid.
	Emergency_Exit1_count = 1
	Emergency_Exit2_count = 1
	while (True):
		C=np.zeros([xyz]) 
		for i in range (Nxyz*cxyz):
			for j in range (Nxyz*cxyz):
				for k in range (Nxyz*cxyz):
					if (RVE_Type=='Gyroid' or RVE_Type=='IWP' or RVE_Type== 'Schwarz' or RVE_Type== 'Diamond' or RVE_Type== 'Neovius' or 
		 					RVE_Type== 'FRD' or RVE_Type== 'FKS' or RVE_Type== 'Octo' or RVE_Type == 'CD', RVE_Type == 'F'):
						
						x=((i*c1)+c1/2) # x-coordinate of center of the voxel     (cxyz*Nxyz)/2.)
						y=((j*c2)+c2/2) # y-coordinate of center of the voxel    -(cxyz*Nxyz)/2.
						z=((k*c3)+c3/2) # z-coordinate of center of the voxel 
					else:
						print("RVE Type is not valid! select: IWP, Schwarz, Gyroid, Diamond, Neovius, FKS or FRD")
						raise AssertionError
					#-------------------------------------------------------------------------------------
					# Calculate F(x,y,z)
					if  (RVE_Type== 'Diamond'):
						a=sin(x)*sin(y)*sin(z)+sin(x)*cos(y)*cos(z)+cos(x)*sin(y)*cos(z)+cos(x)*cos(y)*sin(z) #Diamond
					elif(RVE_Type== 'Gyroid'):
						a=cos(x)*sin(y)+cos(y)*sin(z)+cos(z)*sin(x)  # Gyroide 
					elif (RVE_Type== 'Schwarz'):
						a=cos(x)+cos(y)+cos(z)  #schwarz
					elif (RVE_Type== 'IWP'):  
						a=(cos(2*x)+cos(2*y)+cos(2*z))-2*(cos(x)*cos(y)+cos(y)*cos(z)+cos(z)*cos(x)) #IWP
					elif (RVE_Type== 'Octo'):
						a= 0.6*((cos(x)*cos(y)) + (cos(y)*cos(z)) + (cos(z)*cos(x)))- 0.4*(cos(x)+cos(y)+cos(z)) + 0.25
					elif (RVE_Type== 'F'):
						a= cos(x) * cos(y) * cos(z)
					elif (RVE_Type== 'Neovius'):  
						a=3*(cos(x)+cos(y)+cos(z)) + 4*(cos(x)*cos(y)*cos(z)) #Neovius
					elif (RVE_Type== 'FRD'):  
						a=4*(cos(x)*cos(y)*cos(z)) - (cos(2*x)*cos(2*y)+cos(2*y)*cos(2*z)+cos(2*z)*cos(2*x)) #FRD
					elif (RVE_Type== 'CD'):  
						a = (cos(3*x + y)*cos(z) - sin(3*x -y)*sin(z) + 
                                  cos(x+ 3*y)*cos(z) + sin(x-3*y)*sin(z) + cos(x-y)*cos(3*z) - sin(x+y)*sin(3*z))
					elif (RVE_Type== 'K'):
						a=0.3*(cos(x) + cos(y) + cos(z)) + 0.3*(cos(x)*cos(y) + cos(y)*cos(z)+cos(z)*cos(x)) - 0.4*(cos(2*x) + cos(2*y) + cos(2*z)) + 0.2 
					else:
						a=cos(2*x)*sin(y)*cos(z) + cos(x)*cos(2*y)*sin(z) + sin(x)*cos(y)*cos(2*z) #FKS
					#-------------------------------------------------------------------------------------
					# assign solid phase if F(x,y,z)<t
					if a<=t :
						C[i*cxyz*Nxyz*cxyz*Nxyz+j*cxyz*Nxyz+k]=1
		# Calculate Volume Fraction. 
		VF=100.0*np.sum(C)/len(C)
		# If the volume fraction differs from the target volume fraction by more than 0.1,
		# the t parameter changes accordingly to bring the volume fraction closer to the target volume fraction.
		if (abs(VF-Vol)<.05 or ( (Emergency_Exit1_count>1) and (Emergency_Exit2_count>1) ) ):
			break
		elif (VF<Vol):
			t=t+0.005           #It depends on volCalc.py file (delta in iso_value linspace function)
			Emergency_Exit1_count = Emergency_Exit1_count + 1
		else:
			t=t-0.005
			Emergency_Exit2_count = Emergency_Exit2_count + 1
		#print(VF, t)
	return C, VF, t

#==========================================================
#          		    Part Definition
#==========================================================	
def Part (Nxyz,cxyz,dimscl, ModelName):
	
	NCdmscl = (Nxyz*cxyz)/dimscl
	mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)
	s = mdb.models[ModelName].ConstrainedSketch(name='__profile__', 
		sheetSize=200.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.rectangle(point1=(0.0, 0.0), point2=(NCdmscl, NCdmscl))
	p = mdb.models[ModelName].Part(name='Part-1', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models[ModelName].parts['Part-1']
	p.BaseSolidExtrude(sketch=s, depth=NCdmscl)
	s.unsetPrimaryObject()
	p = mdb.models[ModelName].parts['Part-1']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	del mdb.models[ModelName].sketches['__profile__']
	
	return NCdmscl
#==========================================================
#        				Material
#==========================================================
def Mat (MatParams,gent, ModelName,MatName):

	vacuum_permit = 8.85E-3          #Vacuum Permittivity 
	
	MatNameEndo = MatName[1]
	MatNameExo = MatName[0]

	mdb.models[ModelName].Material(name=MatNameEndo)            # Airy material (Void)
	
	if gent:
	
		endo = mdb.models[ModelName].materials[MatNameEndo]
		endo.Hyperelastic(
			materialType=ISOTROPIC, testData=OFF, type=USER, 
			moduliTimeScale=INSTANTANEOUS, properties=3, table=((MatParams[0], MatParams[1], MatParams[2]), ))
		endo.Conductivity(table=((MatParams[3]*vacuum_permit, ), ))
		mdb.models[ModelName].Material(name=MatNameExo)       # Solid Phase
		exo = mdb.models[ModelName].materials[MatNameExo]
		endo.Hyperelastic(materialType=ISOTROPIC, testData=OFF, type=USER, 
			moduliTimeScale=INSTANTANEOUS, properties=3, table=((MatParams[4], MatParams[5], MatParams[6]), ))
		exo.Conductivity(table=((MatParams[7]*vacuum_permit, ), ))
	
	else:
	
		mdb.models[ModelName].Material(name=MatNameExo)       # Solid Phase
		exo = mdb.models[ModelName].materials[MatNameExo]
		exo.Hyperelastic(materialType=
			ISOTROPIC, table=((MatParams[1], MatParams[2]), ), testData=OFF, type=NEO_HOOKE, 
			volumetricResponse=VOLUMETRIC_DATA)
		exo.Conductivity(table=((MatParams[0]*vacuum_permit, ), ))
		
		endo = mdb.models[ModelName].materials[MatNameEndo]
		endo.Conductivity(table=((MatParams[3]*vacuum_permit, ), 	))
		endo.Hyperelastic(materialType=
			ISOTROPIC, table=((MatParams[4], MatParams[5]), ), testData=OFF, type=NEO_HOOKE, 
			volumetricResponse=VOLUMETRIC_DATA)
	
	mdb.models[ModelName].HomogeneousSolidSection(name='Section-endo', material=MatNameEndo, thickness=None)
	mdb.models[ModelName].HomogeneousSolidSection(name='Section-exo', material=MatNameExo, thickness=None)

	#Assembly--------------------------------------------------------------------

	a = mdb.models[ModelName].rootAssembly
	session.viewports['Viewport: 1'].setValues(displayedObject=a)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(
		optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
	a = mdb.models[ModelName].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models[ModelName].parts['Part-1']
	a.Instance(name='Part-1-1', part=p, dependent=ON)

#==========================================================
#    				Adding mesh, etc. 
#==========================================================
	
def Model (Nxyz,NCdmscl,dimscl, RVE, ModelName,gauss_points):
	
	NC = NCdmscl*dimscl

	if gauss_points == 8:
		elemType1 = mesh.ElemType(elemCode=C3D8HT, elemLibrary=STANDARD)		 # Normal
	else:
		elemType1 = mesh.ElemType(elemCode=C3D8RT, elemLibrary=STANDARD,kinematicSplit=AVERAGE_STRAIN)# 
    		# hourglassControl=ENHANCED)     # Reduced Integration, Hybrid Formulation, Hourglass coontrol 

	a = mdb.models[ModelName].rootAssembly
	p = mdb.models[ModelName].parts['Part-1']
	e = p.edges
	# find edges of the cube
	e1  = e.findAt(((NCdmscl, 0,  .0001*Nxyz),),)
	e2  = e.findAt(((0,  NCdmscl, .0001*Nxyz),),)
	e3  = e.findAt(((0,  0,  .0001*Nxyz),),)
	e4  = e.findAt(((NCdmscl, NCdmscl, .0001*Nxyz),),)
	e5  = e.findAt(((NCdmscl, .0001*Nxyz, 0),),)
	e6  = e.findAt(((0,  .0001*Nxyz, NCdmscl),),)
	e7  = e.findAt(((NCdmscl, .0001*Nxyz, NCdmscl),),)
	e8  = e.findAt(((0,  .0001*Nxyz, 0),),)
	e9  = e.findAt(((.0001*Nxyz,  0, .0),),)
	e10 = e.findAt(((.0001*Nxyz, NCdmscl,  NCdmscl),),)
	e11 = e.findAt(((.0001*Nxyz,  0, NCdmscl),),)
	e12 = e.findAt(((.0001*Nxyz, NCdmscl,  0),),)
	# seed edges
	p.seedEdgeByNumber(edges=(e9+e10+e11+e12), number=int(NC), constraint=FINER)
	p.seedEdgeByNumber(edges=(e5+e6+e7+e8), number=int(NC), constraint=FINER)
	p.seedEdgeByNumber(edges=(e1+e2+e3+e4), number=int(NC), constraint=FINER)
	p.generateMesh()
	# set element type
	c = p.cells
	cells = c.findAt(((.0001*Nxyz, .0001*Nxyz, .0001*Nxyz),),)
	p.setElementType(regions=(cells, ), elemTypes=(elemType1, ))
	a.regenerate()

	#---------------------------------------------
	# Creat two sets for two diffrent materials [0(void) and 1(solid)]
	
	p.PartFromMesh(name='Part-1-mesh-1', copySets=True)
	p1 = mdb.models[ModelName].parts['Part-1-mesh-1']
	del mdb.models[ModelName].parts['Part-1']
	mdb.models[ModelName].parts.changeKey(fromName='Part-1-mesh-1', toName='Part-1')
	a = mdb.models[ModelName].rootAssembly
	del a.features['Part-1-1']
	p = mdb.models[ModelName].parts['Part-1']
	a.Instance(name='Part-1-1', part=p, dependent=ON)
	
	setMat0=[]
	Mat1_nodes = []
	setMat1=[]
	el=p.elements
	
	for i in range (len(RVE)):
		if RVE[i]==0:
			setMat0.append(el[i].label)
		else:
			setMat1.append(el[i].label)
			nodes_el = list(el[i].connectivity)
			for item in nodes_el:
				Mat1_nodes.append(item+1)
	
	Mat1_nodes_unique = list(set(Mat1_nodes))
	nodes = p.nodes
	elset0 = el.sequenceFromLabels(setMat0)
	elset1 = el.sequenceFromLabels(setMat1)
	region1 = p.Set(elements=elset0, name='SetEndo')
	region2 = p.Set(elements=elset1, name='SetExo')
	#Useful_elems = nodes.sequenceFromLabels(Mat1_nodes_unique)
	
	# Assign Sections
	
	p.SectionAssignment(region=region1, sectionName='Section-endo', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
	p.SectionAssignment(region=region2, sectionName='Section-exo', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
	
	model=mdb.models[ModelName]
	del mdb.models['Model-1']
	
	VF_model = (len(setMat1)/ (len(setMat0)+len(setMat1)))*100
	print('VF_Exo (estimated based on mesh): {0} %: '.format(round(VF_model,2)))

	return model,p,nodes,el,elset0,elset1

#==========================================================
#         		    Step Definition
#==========================================================

def Step(model,StepParams,Pre_strch):
			   
	if Pre_strch:
	
		model.CoupledTempDisplacementStep(amplitude=RAMP, cetol=None, 
			creepIntegration=None, deltmx=None, minInc=0.01, initialInc=0.05, maxInc=0.25, 
			maxNumInc=200, name='Step-1', nlgeom=ON, previous='Initial', response=
			STEADY_STATE)
		
		model.CoupledTempDisplacementStep(name='Step-2', 
			previous='Step-1', response=STEADY_STATE, initialInc=StepParams[0],maxNumInc=StepParams[3], deltmx=None, 
			cetol=None,minInc=StepParams[1],maxInc=StepParams[2],
			  creepIntegration=None, amplitude=RAMP)
	else:

		model.CoupledTempDisplacementStep(name='Step-1', 
    		previous='Initial', response=STEADY_STATE, maxNumInc=StepParams[3], 
    		initialInc=StepParams[0], minInc=StepParams[1], maxInc=StepParams[2], deltmx=None, cetol=None, 
    		creepIntegration=None, amplitude=RAMP, nlgeom=ON)
	
	model.fieldOutputRequests['F-Output-1'].setValues(variables=(
		'S', 'E', 'U','LE', 'RF', 'NT','HFL','CF', 'RFLE', 
		'RFL'))
    
#==========================================================
#                Defining Sitable Sets
#==========================================================

def Sets(model,p,nodes,elements,Nxyz,cxyz,dimscl,tol,loding):
	
	NC = cxyz*Nxyz
	NCdmscl = NC/dimscl
	meshSize = 1.0/dimscl
	a = model.rootAssembly

	nodes_BCX =  nodes.getByBoundingBox(xMin =-tol,xMax=tol)
	nodes_BCY =  nodes.getByBoundingBox(yMin=-tol,yMax=tol)
	
	p.Set('ELEMS', elements=elements)
	p.Set('BC_X', nodes=nodes_BCX)
	p.Set('BC_Y', nodes=nodes_BCY)
	
	if loding:
		nodes_BCZ1 =  nodes.getByBoundingBox(zMin=-tol,zMax=tol)
		nodes_BCZ2 =  nodes.getByBoundingBox(zMin = NCdmscl-tol, zMax = NCdmscl+tol)
		p.Set('BC_Z', nodes=(nodes_BCZ1,nodes_BCZ2))

	else:
		nodes_BCZ =  nodes.getByBoundingBox(zMin=-tol,zMax=tol)
		p.Set('BC_Z', nodes=nodes_BCZ)
	
	load_X = elements.getByBoundingBox(xMin = NCdmscl-(3*meshSize/2), xMax = NCdmscl+(3*meshSize/2))
	#PreS_Z = elements.getByBoundingBox(zMin = NCdmscl-(3*meshSize/2), zMax = NCdmscl+(3*meshSize/2))
		
	#nodes_PreS_X = elems.getByBoundingBox(xMin = CN-(3*meshSize/2), xMax = CN+(3*meshSize/2))
	#elems_PreS_Z = elems.getByBoundingBox(zMin = CN-(3*meshSize/2), zMax = CN+(3*meshSize/2))
	
	#p.Set('PreS_X', elements=PreS_X)
	#p.Set('PreS_Z', elements=PreS_Z)

	nodes_U1 =  nodes.getByBoundingBox(xMin =NCdmscl-tol,xMax=NCdmscl+tol)
	nodes_U2 =  nodes.getByBoundingBox(yMin=NCdmscl-tol,yMax=NCdmscl+tol)
	nodes_U3 =  nodes.getByBoundingBox(zMin=NCdmscl-tol,zMax=NCdmscl+tol)
	
	#2. In order to track avg Reaction Forces(RF) & Displacements(U) on surfaces, the following sets are defined.
	
	p.Set('U1', nodes=nodes_U1)
	p.Set('U2', nodes=nodes_U2)
	p.Set('U3', nodes=nodes_U3)

	HFregion = elements.getByBoundingBox(yMin =NCdmscl-meshSize-tol,yMax=NCdmscl+meshSize+tol)
	p.Surface(face6Elements=HFregion, name='Surf-HF')
	p.Surface(face2Elements=load_X, name='Load-X')
	
	a.regenerate()

#========================================================================
#          				Boundary Conditions
#========================================================================

def BCs(model,Vlt_Bot,HFL_mgn,Pre_strch,fLoad,ModelName):   

	a = mdb.models[ModelName].rootAssembly

	region_BCX = a.instances['Part-1-1'].sets['BC_X']
	region_BCY = a.instances['Part-1-1'].sets['BC_Y']
	region_BCZ = a.instances['Part-1-1'].sets['BC_Z']
	region_loadX = a.instances['Part-1-1'].surfaces['Load-X']
	BC_U1 = a.instances['Part-1-1'].sets['U1']
	BC_U3 = a.instances['Part-1-1'].sets['U3']
	SURF_HF = a.instances['Part-1-1'].surfaces['Surf-HF']
	
	#========================================================================
	
	model.DisplacementBC(name='X-Constraint', createStepName='Initial', 
	region=region_BCX, u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	model.DisplacementBC(name='Y-Constraint', createStepName='Initial', 
	region=region_BCY, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	model.DisplacementBC(name='Z-Constraint', createStepName='Initial', 
	region=region_BCZ, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
	
	#========================================================================

	model.TemperatureBC(name='Voltage-Bottom', createStepName='Initial', 
	region=region_BCY, distributionType=UNIFORM, fieldName='', magnitude=Vlt_Bot)    

	#model.TemperatureBC(name='Voltage-Up', createStepName='Step-1', 
	#	region=region_Y_Up, distributionType=UNIFORM, fieldName='', magnitude=Vlt_Up)   

	if Pre_strch:         
		
		model.Pressure(name='Load_X', createStepName='Step-1', 
			region=region_loadX, distributionType=UNIFORM, field='', magnitude=-fLoad, 
			amplitude=UNSET)
		
		model.DisplacementBC(name='Z-Constraint2', createStepName='Initial', 
			region=BC_U3, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
			amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
		
		model.SurfaceHeatFlux(name='SHF', createStepName='Step-2', 
    		region=SURF_HF, magnitude=HFL_mgn)
	else:
		model.SurfaceHeatFlux(name='SHF', createStepName='Step-1', 
		    region=SURF_HF, magnitude=HFL_mgn)

#====================================================================
#       				Job Definition
#====================================================================

def Job (tempInpName,modelName):
	assem = mdb.models[modelName].rootAssembly
	assem.regenerate()
	mdb.Job(name=tempInpName, model=modelName)
	assem.regenerate()
	mdb.jobs[tempInpName].writeInput(consistencyChecking=OFF)

#====================================================================
#              				Details
#====================================================================

def Detail(MatParams,SimParams,nu,VF,Pre_strch,fLoad,Nxyz,cxyz,dimScl,ModelName,MatName):    #region_X,region_Y,region_Z)
	
	a = mdb.models[ModelName].rootAssembly
	a.regenerate()
	rootAssem_ins = a.instances['Part-1-1']
	num_nodes = len(rootAssem_ins.nodes)                 # Number of nodes is stored in this container and returned through function
	num_elems = len(rootAssem_ins.elements)              # Number of elements is stored in this container and returned through function
	NCDimscl = Nxyz*cxyz/dimScl
	vac_permit  = 8.85E-3
	
	with open(ModelName + '_dtl.txt', 'a') as f:          # creating a csv writer object
		f.write('-------------------------------\n')
		f.write('-------- Abaqus Detail --------\n')
		f.write('-------------------------------\n\n')
		f.write('Model Name:\t' + ModelName + '\n' )
		f.write('Type of elements would be Quadilaterial \n')
		f.write('Num of nodes would be: {0}\n'.format(num_nodes))
		f.write('Num of elements would be: {0}\n'.format(num_elems))
		f.write('Nominal Volume Fraction of Exo phase: {0}%\n'.format(VF[0]))
		f.write('Actual Volume Fraction of Exo phase: {0}%\n'.format(round(VF[1],2)))
		f.write('Num of element gaussian points: {0}\n'.format(SimParams[2]))
		f.write('\n--- %% Geometry Parameters %% ---\n\n')
		f.write('Num of Cells (each dir): {0}\n'.format(Nxyz))
		f.write('Num of tesselation (each dir):  {0}\n'.format(cxyz))
		f.write('Scale Factor: {0} \n'.format(dimScl))
		f.write('Length of RVE (each dir): {0}  mm\n'.format(NCDimscl))

		if Pre_strch:
			f.write('Lateral Pre-Stress: {0} kPa\n'.format(fLoad))
		else:
			f.write('Lateral Pre-Stretch : {0} (without Pre-stretch)\n'.format(1))
		
		f.write('\n--- %% Material Parameters %% ---\n\n')
		f.write('Shear Modulus of Exo-Skeleton:  {0}\tkPa\n'.format(MatParams[1]*2))
		f.write('Exo-Skeleton material:  {0}\t\n'.format(MatName[0])) 
		f.write('Relative Permittiviy of Exo-Skeleton:  {0}\n'.format(MatParams[0]))
		f.write('Permittiviy of Exo-Skeleton:  {0}\n'.format(vac_permit*MatParams[0]))
		f.write('Endo-Skeleton material:  {0}\t\n'.format(MatName[1]))   
		f.write('Shear Modulus of Endo-Skeleton:  {0}\tkPa\n'.format(MatParams[4]*2)) 
		f.write('Relative Permittiviy of Endo-Skeleton:  {0}\n'.format(MatParams[3]))
		f.write('Permittiviy of Exo-Skeleton:  {0}\n'.format(vac_permit*MatParams[3]))
		f.write('Poisson ratio of both phases:  {0}\n'.format(nu))

		f.write('\n--- %% Self-consistent Unit %% ---\n\nLength:   mm\nForce:   mN\n'+
			'Voltage: kV\nPermittivity: pF/mm\nStress(Pressure):  kPa\n'+
			'Time:  s\nMass:  Kg\nDensity:  Kg/mm3\n\n')
	
#=================================================================================

def flatten_extend(matrix):
	flat_list = []
	for row in matrix:
		flat_list.extend(row)
	return flat_list

#=================================================================================

def phaseLists(InpFile):

	findPhrase1 = '*Elset, elset=SetEndo'
	findPhrase2 = '*Elset, elset=SetExo'

	toggle1 = 0 		#=the first phrase is found
	#toggle2 = 0 		#=the second phrase is not found

	tempInp = open(InpFile, 'r')

	#coupledInp = open(Inpfile_exp, 'w')

	while toggle1==0: 	#=continue until the both phrases are found
		line = tempInp.readline()
		#coupledInp.write(line)
		if line.startswith(findPhrase1): toggle1 = 1	#the first phrase is found
		#if toggle1 == 1 and line.startswith(findPhrase2): toggle2 = 1	#the second phrase is found

	VoidElems = []
	ElasElems = []

	line = tempInp.readline()

	while line.startswith(findPhrase2)!=True:	#continue until next command starting with '*'
		elem = re.findall('\d+', line)
		elem = [int(x) for x in elem]
		VoidElems.append(elem)
		line = tempInp.readline()

	while line.startswith('*Nset')!= True:
		elem = re.findall('\d+', line)
		elem = [int(x) for x in elem]
		ElasElems.append(elem)
		line = tempInp.readline()
	
	
	return VoidElems, ElasElems

#==================================================================
# 				        Running the Model
#==================================================================

def Run(MatParams,SimParams,uelFile,Voidelems_flat,RunNames):
	
	vacuum_permit = 8.85E-3            #Vacuum Permittivity 
	
	exo_permit= vacuum_permit * MatParams[0]
	endo_permit= vacuum_permit * MatParams[3]

	line_uelDef = '*user element, nodes = {0}, coordinates = 3, type = U1, \n'.format(SimParams[1])
	line_uelProp_exo = '{0}, {1} \n'.format(exo_permit,SimParams[2])
	line_uelProp_endo = '{0}, {1} \n'.format(endo_permit,SimParams[2])

	uelDefinition1 = ('**\n'+
					'**\n'+
					'**\n'+
					line_uelDef+
					'properties = 1, i properties = 1, variables = 1, unsymm\n'+
					'1, 2, 3, 11\n')
	uelElemsInfo_exo = '*element, type = U1, elset = uels_exo\n'
	uelelemsInfo_endo = '*element, type = U1, elset = uels_endo\n'

	uelProperties_exo = ('*uel property, elset = uels_exo\n'+
					line_uelProp_exo+
					'**\n'+
					'**\n'+
					'**\n')


	uelProperties_endo = ('*uel property, elset = uels_endo\n'+
					line_uelProp_endo+
					'**\n'+
					'**\n'+
					'**\n')

	# --------------------------------------------------------------

	coupledInp = open(RunNames[1]+'.inp', 'w')
	
	coupledInpName = RunNames[1]

	Voidelems_flat_np = np.array(Voidelems_flat)

	tempInp = open(RunNames[0]+'.inp', 'r')

	#tempInp.seek(0, 0)
	
	startUelElemNum_endo = 10000000
	startUelElemNum_exo = 20000000

	toggle2 = 0
	findPhrase_elem = '*Element'
	
	while toggle2==0: 	#=continue until the both phrases are found
		line = tempInp.readline()
		#print(line)
		coupledInp.write(line)
		if line.startswith(findPhrase_elem): toggle2 = 1	#the first phrase is found

	elemUel_air=''
	elemUel_elas='' 
	line = tempInp.readline()
	#elemLine = tempInp.readline()
	#print(elemLine)
	while line.startswith('*Elset,')!=True:	#continue until next command starting with '*'
		elem = re.findall('\d+', line)
		elem = [int(x) for x in elem]
		coupledInp.write(line)
		#print(line)
		#print(np.isin(elem[0],Voidelems_flat_np))
		if np.isin(elem[0],Voidelems_flat_np):
			elem[0] = elem[0] + startUelElemNum_endo
			elemUel_air = elemUel_air+',   '.join(map(str,elem))+'\n'
		else:
			elem[0] = elem[0] + startUelElemNum_exo
			elemUel_elas = elemUel_elas+',   '.join(map(str,elem))+'\n'
		
		line = tempInp.readline()
		
		#elem[0] = elem[0] + startUelElemNum
		#elemUel = elemUel+',   '.join(map(str,elem))+'\n'
		#elemLine = tempInp.readline()


	coupledInp.write(uelDefinition1)
	coupledInp.write(uelelemsInfo_endo)
	coupledInp.write(elemUel_air)
	coupledInp.write(uelProperties_endo)

	coupledInp.write(uelElemsInfo_exo)
	coupledInp.write(elemUel_elas)
	coupledInp.write(uelProperties_exo)

	coupledInp.write(line)
	lines = tempInp.read()
	coupledInp.write(lines)
	tempInp.close()
	coupledInp.close()
	
	# Creating a job file
	
	mdb.JobFromInputFile(name=RunNames[1], 
    	inputFileName=RunNames[1]+'.inp',
    	userSubroutine=uelFile)
	mdb.jobs[coupledInpName].setValues(numDomains=SimParams[0], 
    	activateLoadBalancing=False, 
    	numCpus=SimParams[0])
	
	#time1om = time.time()
	mdb.jobs[coupledInpName].submit()
	
	#mdb.JobFromInputFile(name=RunNames[1], 
	#	inputFileName=RunNames[1]+'.inp')
	#mdb.jobs[coupledInpName].setValues(getMemoryFromAnalysis=True, 
	#	explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
	#	userSubroutine=uelFile, 
	#	scratch='', numThreadsPerMpiProcess=1, 
	#	multiprocessingMode=DEFAULT, numCpus=SimParams[0], numDomains=2*SimParams[0], numGPUs=1)

	#mdb.jobs[coupledInpName].waitForCompletion()

	#time2om = time.time()
	#with open(RunNames[1]+'_dtl.txt', 'a') as f:          # creating a csv writer object
	#	f.write('--- %% Run Information %% ---')
	#	f.write('\nElapsed Time (in mins): {0}\n'.format((time2om-time1om)/60.0))