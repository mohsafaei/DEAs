#=================================================================================================
# In silico actuation performance investigation of dielectric elastomers with TPMS geometries
#### Link: https://doi.org/10.1016/j.euromechsol.2024.105540   ####
#=================================================================================================

# Part of the code for creating TPMS structures was developed by Fayyaz Nosouhi (dehnavifn@gmail.com) and Saeed Khaleghi (saeedkhaleghi123@gmail.com).
# The remaining parts of the code, including the definition of DEA composites and their implementation in ABAQUS, were developed by Mohammad Ali Safaei (mohammadsf1998@gmail.com).
# This Python script requires the use of a UEL subroutine, which implements the constitutive equations for a DEA element. In this regard, the UEL subroutine developed by Ehsan Hajiesmaili was utilized. 
# University of Tehran, December 2024

#=================================================================================

# ------------------------------------
#             Modules
# ------------------------------------

import mech_TPMS
from mech_TPMS import MyRVE
from mech_TPMS import Part
from mech_TPMS import Mat
from mech_TPMS import Model
from mech_TPMS import Step
from mech_TPMS import Sets
from mech_TPMS import BCs
from mech_TPMS import Job
from mech_TPMS import flatten_extend
from mech_TPMS import Detail
from mech_TPMS import Run
from mech_TPMS import phaseLists
import shutil
import os
import time
from operator import itemgetter 
import numpy as np

# --------------------------------------
#         Geometry Params
# --------------------------------------

# Num. of Voxels in each direction
Nxyz = 30    # Nx = Ny = Nz

# Number of cells in each direction
cxyz = 4   # cx = cy = cz 

# Dimension Scaling in order to bring voltages lower
Dim_Scl = cxyz*Nxyz*10

l_RVE = Nxyz*cxyz/Dim_Scl

# Volume Fraction of Exo Phase (Solid)
VF_nom = 50 # Volume Fraction (%) [0,100]

# RVE_Type
RVE_Type = 'Schwarz'  # [IWP,   Schwarz,  Gyroid,  Diamond, Neovius, FKS, FRD, Octo,CD]

Pre_strch = False
fLoad = 150   #kPa    60.0

PreS_Coff = 1.25
PreS_Value = (PreS_Coff-1.00)*l_RVE

vlt_bot = 0.0
HFL_mgn  = 10.0



ModelName= RVE_Type+'_'+str(VF_nom)+'_'+str(Nxyz)+'_'+str(Dim_Scl) #+'_f150'

ModelName=ModelName.replace('.', '_')

# ------------------------------------
#         Material Parameters
# ------------------------------------

Gent = False
mu_exo = 100.0          #kPa  Shear modulus of Eco-Flex   
mu_endo = 440.0          #kPa  Shear modulus of Sylgard 10_1
mu_comb = (mu_exo*VF_nom/100) + ((100-VF_nom)/100*mu_endo)

nu = 0.475          #     Poisson ratio for nearly incompressible material   
                    #     [0.452,0.475,0.490,0.495,0.4995,0.49995]

c10_exo = mu_exo/2          #kPa                   ABAQUS Parameters for Neo-Hookean material      
c10_endo = mu_endo/2          #kPa           
D1_exo  = (3 *(1- 2 * nu)/(2*c10_exo*(1+nu)) )     #1/kPa     ABAQUS Parameters ....
D1_endo  = (3 *(1- 2 * nu)/(2*c10_endo*(1+nu)) )     #1/kPa     

rlt_permit_exo = 3.2             #Relative Permittivity of Eco-flex
rlt_permit_endo = 2.7            #Relative Permittivity of Sylgard 10_1
rlt_permit_comb = (rlt_permit_exo*VF_nom/100) + ((100-VF_nom)/100*rlt_permit_endo)

MatParams = [rlt_permit_exo,c10_exo,D1_exo,rlt_permit_endo,c10_endo,D1_endo,rlt_permit_comb,mu_comb]

MatName = ['Sylgard 10_1','Eco-flex']
# ------------------------------------
#         Step Parameters
# ------------------------------------

InitInc = 0.001
minInc = 1E-04
maxInc = 0.005
maxNumInc = 100000
StepParams = [InitInc,minInc,maxInc,maxNumInc]

# --------------------------------------------
#         Simulation Parameters
# --------------------------------------------

uelFile = 'Static.for'
ncpu = int(input('Nubmber of CPUs for this job: '))
Num_nodes = 8
gauss_points = 1

SimParams = [ncpu,Num_nodes,gauss_points]
tempInpName = 'temp_' + ModelName       #Creating a temp inp file to read it and rewrite it simultaneously
coupledINPName =  ModelName

RunNames = [tempInpName,coupledINPName,ModelName]

# -------------------------------------------------------------
tol = 0.001
RVE,VF_act,LSP_t = MyRVE(Nxyz, cxyz, RVE_Type, VF_nom)
print('Level Set Parameter(t): {0}\n & VF_actual: {1} %:'.format(LSP_t,round(VF_act,2)))

# Creating FEM Model
VF = [VF_nom,VF_act]

NCdmscl = Part(Nxyz,cxyz,Dim_Scl,ModelName)
Mat(MatParams,Gent,ModelName,MatName)
model,p,nodes,elements,elset0,elset1 = Model(Nxyz,NCdmscl,Dim_Scl, RVE, ModelName,gauss_points)
Step(model,StepParams,Pre_strch)
Sets(model,p,nodes,elements,Nxyz,cxyz,Dim_Scl,tol,Pre_strch)
BCs(model,0.0,HFL_mgn,Pre_strch,fLoad,ModelName)
Job(tempInpName,ModelName)
VoidElems,ElasElems = phaseLists(tempInpName+'.inp')
Voidelems_flat = flatten_extend(VoidElems)
Elaselems_flat = flatten_extend(ElasElems)
Detail(MatParams,SimParams,nu,VF,Pre_strch,fLoad,Nxyz,cxyz,Dim_Scl,ModelName,MatName)
Run(MatParams,SimParams,uelFile,Voidelems_flat,RunNames)
