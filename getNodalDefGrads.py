#Property of Not Real Engineering
#Copyright 2020 Not Real Engineering - All Rights Reserved You may not use,
#           distribute and modify this code without the written permission
#           from Not Real Engineering.
############################################################################
##             Reading the ODB file                                       ##
############################################################################


from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
from odbAccess import openOdb
import odbAccess
import math
import numpy
import os        # Operating system
import shutil    # copying or moving files

#Open text file to write results
# root_dir = 'C:/Users/aj686/OneDrive - University Of Cambridge/Documents/'


work_dir = 'C:/Users/aj686/OneDrive - University Of Cambridge/Documents/HexElverify'
os.chdir(work_dir)
odbname='Job-5-1'         # set odb name here
path=work_dir+'/'                    # set odb path here (if in working dir no need to change!)
myodbpath=path+odbname+'.odb'
odb=openOdb(myodbpath)

tol = 1e-6
toly = 1e-5

set = odb.rootAssembly.instances['PART-1'].nodeSets['SET-1']
elset = odb.rootAssembly.instances['PART-1'].elementSets['SET-1']
All_disp = odb.steps['Step-1'].frames[1].fieldOutputs['U'].values

x1 = []
x2 = []
x3 = []
x4 = []

u1 = []
u2 = []
u3 = []

connectivity = numpy.zeros((len(elset.elements),len(elset.elements[0].connectivity)),dtype=numpy.int64)

# sortie = open(path+odbname+'_NodPosDisps.csv' , 'w')
# sortie.write('x1,x2,x3,u1,u2,u3\n')

for q in range(len(set.nodes)):
	x1.append(set.nodes[q].coordinates[0])
	x2.append(set.nodes[q].coordinates[1])
	x3.append(set.nodes[q].coordinates[2])
	x4.append(set.nodes[q].label)
	u1.append(All_disp[q].data[0])
	u2.append(All_disp[q].data[1])
	u3.append(All_disp[q].data[2])
	# sortie.write('%f,%f,%f,%f,%f,%f\n'%(set.nodes[q].coordinates[0],set.nodes[q].coordinates[1],set.nodes[q].coordinates[2],All_disp[q].data[0],All_disp[q].data[1],All_disp[q].data[2]))

# sortie.close()
#
# sortie = open(path+odbname+'_elems.csv' , 'w')
# sortie.write('node1,node2,node3,node4\n')

for q in range(len(elset.elements)):
	for k in range(len(elset.elements[0].connectivity)):
		connectivity[q,k] = (elset.elements[q].connectivity[k])
# 		sortie.write('%d'%(connectivity[q,k]))
# 		if k<len(elset.elements[0].connectivity)-1:
# 			sortie.write(',')
# 	sortie.write('\n')
# sortie.close()

x1 = numpy.array(x1)
x2 = numpy.array(x2)
x3 = numpy.array(x3)
x4 = numpy.array(x4)

u1 = numpy.array(u1)
u2 = numpy.array(u2)
u3 = numpy.array(u3)

# u1[169] =numpy.nan
# u2[169] =numpy.nan
# u3[169] =numpy.nan

def getGradTet(nodeslist):
	p1 = nodeslist[0]
	p2 = nodeslist[1]
	p3 = nodeslist[2]
	p4 = nodeslist[3]
	ppVol = numpy.dot(numpy.cross(p2-p1,p3-p1),(p4-p1))
	return numpy.row_stack((numpy.cross(p4-p2,p3-p2)/ppVol,numpy.cross(p3-p1,p4-p1)/ppVol,numpy.cross(p4-p1,p2-p1)/ppVol,numpy.cross(p2-p1,p3-p1)/ppVol)) , ppVol

def getGradHex(nodeslist):
	p1 = nodeslist[0]
	p2 = nodeslist[1]
	p3 = nodeslist[2]
	p4 = nodeslist[3]
	p5 = nodeslist[4]
	p6 = nodeslist[5]
	p7 = nodeslist[6]
	p8 = nodeslist[7]

	xi = 0
	eta = 0
	beta = 0

	X = numpy.column_stack((p1,p2,p3,p4,p5,p6,p7,p8))
	N = 1/8*numpy.array([(1-xi)*(1+eta)*(1+beta),
	(1-xi)*(1-eta)*(1+beta),
	(1+xi)*(1-eta)*(1+beta),
	(1+xi)*(1+eta)*(1+beta),
	(1-xi)*(1+eta)*(1-beta),
	(1-xi)*(1-eta)*(1-beta),
	(1+xi)*(1-eta)*(1-beta),
	(1+xi)*(1+eta)*(1-beta)])
	dNdxi = 1/8*numpy.array([[-(1+eta)*(1+beta), (1-xi)*(1+beta), (1-xi)*(1+eta)],
	[-(1-eta)*(1+beta), -(1-xi)*(1+beta), (1-xi)*(1-eta)],
	[(1-eta)*(1+beta), -(1+xi)*(1+beta), (1+xi)*(1-eta)],
	[(1+eta)*(1+beta), (1+xi)*(1+beta), (1+xi)*(1+eta)],
	[-(1+eta)*(1-beta), (1-xi)*(1-beta), -(1-xi)*(1+eta)],
	[-(1-eta)*(1-beta), -(1-xi)*(1-beta), -(1-xi)*(1-eta)],
	[(1-eta)*(1-beta), -(1+xi)*(1-beta), -(1+xi)*(1-eta)],
	[(1+eta)*(1-beta), (1+xi)*(1-beta), -(1+xi)*(1+eta)]])

	J = numpy.matmul(X,dNdxi)
	detJ = numpy.linalg.det(J)
	ppVol = 8*detJ
	Jinv = numpy.linalg.inv(J)
	dNdx = nump.matmul(dNdxi,Jinv)
	return dNdx , ppVol


sortie = open(path+odbname+'_DefGrads.csv' , 'w')
sortie.write('F11,F22,F33,F12,F13,F23,F21,F31,F32,elVol\n')#Just to conform to Abaqus' silly format

for q in range(len(elset.elements)):
	nodeslist = []
	ulist = []
	for count in range(len(elset.elements[0].connectivity)):
		seek= (x4==connectivity[q,count])
		nodeslist.append(numpy.array([x1[seek],x2[seek],x3[seek]]).squeeze())
		ulist.append(numpy.array([u1[seek],u2[seek],u3[seek]]).squeeze())

	p1 = nodeslist[0]
	p2 = nodeslist[1]
	p3 = nodeslist[2]
	p4 = nodeslist[3]
	p5 = nodeslist[4]
	p6 = nodeslist[5]
	p7 = nodeslist[6]
	p8 = nodeslist[7]

	xi = 0
	eta = 0
	beta = 0

	X = numpy.column_stack((p1,p2,p3,p4,p5,p6,p7,p8))
	N = 1/8*numpy.array([(1-xi)*(1+eta)*(1+beta),
	(1-xi)*(1-eta)*(1+beta),
	(1+xi)*(1-eta)*(1+beta),
	(1+xi)*(1+eta)*(1+beta),
	(1-xi)*(1+eta)*(1-beta),
	(1-xi)*(1-eta)*(1-beta),
	(1+xi)*(1-eta)*(1-beta),
	(1+xi)*(1+eta)*(1-beta)])
	dNdxi = numpy.array([[-(1+eta)*(1+beta), (1-xi)*(1+beta), (1-xi)*(1+eta)],
	[-(1-eta)*(1+beta), -(1-xi)*(1+beta), (1-xi)*(1-eta)],
	[(1-eta)*(1+beta), -(1+xi)*(1+beta), (1+xi)*(1-eta)],
	[(1+eta)*(1+beta), (1+xi)*(1+beta), (1+xi)*(1+eta)],
	[-(1+eta)*(1-beta), (1-xi)*(1-beta), -(1-xi)*(1+eta)],
	[-(1-eta)*(1-beta), -(1-xi)*(1-beta), -(1-xi)*(1-eta)],
	[(1-eta)*(1-beta), -(1+xi)*(1-beta), -(1+xi)*(1-eta)],
	[(1+eta)*(1-beta), (1+xi)*(1-beta), -(1+xi)*(1+eta)]],dtype=numpy.double)*0.125

	J = numpy.matmul(X,dNdxi)
	detJ = numpy.linalg.det(J)
	ppVol = 8*detJ
	# breakpoint()
	Jinv = numpy.linalg.inv(J)
	integrator = numpy.matmul(dNdxi,Jinv)
	# integrator,ppVol=getGradHex(nodeslist) #provides a row of integrator table
	ulist = numpy.array(ulist)
	ulist = ulist.transpose()
	F = numpy.matmul(ulist,integrator) + numpy.eye(3)
	sortie.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n'%(F[0,0],F[1,1],F[2,2],F[0,1],F[0,2],F[1,2],F[1,0],F[2,0],F[2,1],ppVol))

sortie.close()
	#Property of Not Real Engineering
