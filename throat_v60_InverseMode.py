import Sofa.Core
import SofaRuntime
import numpy as np
from math import cos,atan,acos,sqrt,sin,pi
from splib3.numerics import Vec3, Quat,sdiv
from scipy.spatial.transform import Rotation as R
from splib3.objectmodel import setData
from splib3.animation import animate, AnimationManager
import os
import time
#Screenshot
import win32gui
# from PyQt5.QtWidgets import QApplication
# from PyQt5.QtGui import *
import sys
# from ArduinoInterface_730 import Controller


path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'
dirPath = os.path.dirname(os.path.abspath(__file__))+'/'
Img_path = os.path.dirname(os.path.abspath(__file__))+'/img/'
##########################################
# Some Settings                          #
##########################################
#Setting CableActuator position:
#first point radius
length1 = 2.2
#other point radius
cablelength=2.2
#Drawing Circle: 
radius=40
modelX=0
modelY=0
modelZ=0
#TargetPoint position in inverse mode: 
target_position=[0+modelX, -55+modelY,100+modelZ]
# target_position=[0, -55,160]

#EffectorPoint on the robot:
robot_decisionpoint=[[0., 0., 130]]
#Static patient position
patient_position=[0,-50,140,0,0,0,1]
print_flag=0
#For ScreenShot
hwnd_title = dict()

##########################################
# CableGeometry positions of every point #
##########################################
cableGeometry1 = [[ 0,-length1, 0.],[0,-cablelength, 2.2],[0,-cablelength, 3],[0,-cablelength, 4],[0,-cablelength, 6],[0,-cablelength,10],[0,-cablelength, 15],
                  [0,-cablelength, 20],[0,-cablelength, 25],[0,-cablelength,30],[0,-cablelength,35],[0,-cablelength, 40],[0,-cablelength,45],
                  [0,-cablelength,50],[0,-cablelength,55],[0,-cablelength, 58],[0,-cablelength, 60]]
cableGeometry2 = [[ -length1*cos(pi/6),length1*cos(pi/3), 0.],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 2.2],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 3],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 4],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 6],[-cablelength*cos(pi/6),cablelength*cos(pi/3),10],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 15],
                  [-cablelength*cos(pi/6),cablelength*cos(pi/3), 20],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 25],[-cablelength*cos(pi/6),cablelength*cos(pi/3),30],[-cablelength*cos(pi/6),cablelength*cos(pi/3),35],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 40],[-cablelength*cos(pi/6),cablelength*cos(pi/3),45],
                  [-cablelength*cos(pi/6),cablelength*cos(pi/3),50],[-cablelength*cos(pi/6),cablelength*cos(pi/3),55],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 58],[-cablelength*cos(pi/6),cablelength*cos(pi/3), 60]]
cableGeometry3 = [[ length1*cos(pi/6),length1*cos(pi/3), 0.],[cablelength*cos(pi/6),cablelength*cos(pi/3), 2.2],[cablelength*cos(pi/6),cablelength*cos(pi/3), 3],[cablelength*cos(pi/6),cablelength*cos(pi/3), 4],[cablelength*cos(pi/6),cablelength*cos(pi/3), 6],[cablelength*cos(pi/6),cablelength*cos(pi/3),10],[cablelength*cos(pi/6),cablelength*cos(pi/3), 15],
                  [cablelength*cos(pi/6),cablelength*cos(pi/3), 20],[cablelength*cos(pi/6),cablelength*cos(pi/3), 25],[cablelength*cos(pi/6),cablelength*cos(pi/3),30],[cablelength*cos(pi/6),cablelength*cos(pi/3),35],[cablelength*cos(pi/6),cablelength*cos(pi/3), 40],[cablelength*cos(pi/6),cablelength*cos(pi/3),45],
                  [cablelength*cos(pi/6),cablelength*cos(pi/3),50],[cablelength*cos(pi/6),cablelength*cos(pi/3),55],[cablelength*cos(pi/6),cablelength*cos(pi/3), 58],[cablelength*cos(pi/6),cablelength*cos(pi/3), 60]]
cableGeometry4 = [[0, length1, 0.],[ 0,cablelength, 2.2],[ 0,cablelength, 3],[ 0,cablelength, 4],[ 0,cablelength,6],[ 0,cablelength,10],[ 0,cablelength, 15],
                  [0,cablelength, 20],[0,cablelength, 25],[0,cablelength,30],[0,cablelength,35],[0,cablelength, 40],[ 0,cablelength,45],
                  [0,cablelength,50],[0,cablelength,55],[0,cablelength, 58],[0,cablelength,60],[ 0,cablelength, 64],[ 0,cablelength, 67],[ 0,cablelength, 70],[ 0,cablelength, 73],[ 0,cablelength,76],[ 0,cablelength, 82],
                  [0,cablelength, 86],[0,cablelength, 90],[0,cablelength,94],[0,cablelength,98],[0,cablelength, 104],[ 0,cablelength,109],
                  [0,cablelength,112],[0,cablelength,116],[0,cablelength, 120],[0,cablelength, 124],[0,cablelength, 127],[ 0, cablelength,129]]
cableGeometry6 = [[length1*cos(pi/6), -length1*cos(pi/3), 0.],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 2.2],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 3],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 4],[ cablelength*cos(pi/6),-cablelength*cos(pi/3),6],[ cablelength*cos(pi/6),-cablelength*cos(pi/3),10],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 15],
                  [cablelength*cos(pi/6),-cablelength*cos(pi/3), 20],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 25],[cablelength*cos(pi/6),-cablelength*cos(pi/3),30],[cablelength*cos(pi/6),-cablelength*cos(pi/3),35],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 40],[ cablelength*cos(pi/6),-cablelength*cos(pi/3),45],
                  [cablelength*cos(pi/6),-cablelength*cos(pi/3),50],[cablelength*cos(pi/6),-cablelength*cos(pi/3),55],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 58],[cablelength*cos(pi/6),-cablelength*cos(pi/3),60],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 64],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 67],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 70],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 73],[cablelength*cos(pi/6),-cablelength*cos(pi/3),76],[ cablelength*cos(pi/6),-cablelength*cos(pi/3), 82],
                  [cablelength*cos(pi/6),-cablelength*cos(pi/3), 86],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 90],[cablelength*cos(pi/6),-cablelength*cos(pi/3),94],[cablelength*cos(pi/6),-cablelength*cos(pi/3),98],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 104],[ cablelength*cos(pi/6),-cablelength*cos(pi/3),109],
                  [cablelength*cos(pi/6),-cablelength*cos(pi/3),112],[cablelength*cos(pi/6),-cablelength*cos(pi/3),116],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 120],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 124],[cablelength*cos(pi/6),-cablelength*cos(pi/3), 127],[ cablelength*cos(pi/6),-cablelength*cos(pi/3),129]]
cableGeometry5= [[-length1*cos(pi/6), -length1*cos(pi/3), 0.],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 2.2],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 3],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 4],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),6],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),10],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 15],
                  [-cablelength*cos(pi/6),-cablelength*cos(pi/3), 20],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 25],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),30],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),35],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 40],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),45],
                  [-cablelength*cos(pi/6),-cablelength*cos(pi/3),50],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),55],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 58],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),60],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 64],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 67],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 70],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 73],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),76],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3), 82],
                  [-cablelength*cos(pi/6),-cablelength*cos(pi/3), 86],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 90],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),94],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),98],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 104],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),109],
                  [-cablelength*cos(pi/6),-cablelength*cos(pi/3),112],[-cablelength*cos(pi/6),-cablelength*cos(pi/3),116],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 120],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 124],[-cablelength*cos(pi/6),-cablelength*cos(pi/3), 127],[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),129]]

##########################################
# Setting a target in Inverse mode       #
##########################################
def effectorTarget(parentNode, position=target_position):
    target = parentNode.addChild('Target')
    target.addObject('EulerImplicitSolver', firstOrder=True)
    target.addObject('CGLinearSolver')
    target.addObject('MechanicalObject', name='dofs', position=position, showObject=True, showObjectScale=2, drawMode=2, showColor=[1., 1., 1., 1.])
    target.addObject('UncoupledConstraintCorrection')
    return target
####################
# Throat                                 #
##########################################
class Throat():
    #Setting basic properties about your throat robot
    def __init__(self, parentNode, youngModulus=10000, poissonRatio=0.45, totalMass=0.07):  
        self.node = parentNode.addChild('Throat')
        #Loading vtk file and some necessary steps
        self.node.addObject('MeshVTKLoader', name='loader', filename=path+'728_v1.vtk')
        self.node.addObject('MeshTopology', src='@loader', name='container')
        self.node.addObject('MechanicalObject', name='dofs', template='Vec3')
        self.node.addObject('UniformMass', totalMass=totalMass)
        self.node.addObject('TetrahedronFEMForceField', template='Vec3', name='FEM', method='large', poissonRatio=poissonRatio,  youngModulus=youngModulus)
        #Adding cables to your robot
        self.__addCables()
        self.addTip()
        self.addFixedbox()

    def __addCables(self):
        #Define a cable      
        #1.add a child node
        cable1 = self.node.addChild("cable1")
        #2.create a mech object and set position
        cable1.addObject("MechanicalObject", name='dof1', position=cableGeometry1)
        #3.add actuator for inverse mode and set params
        cable1.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry1))),
                        maxPositiveDisp='60',
                        maxDispVariation='1.5',
                        hasPullPoint='0',
                        minForce=0)
        #4.mapping
        cable1.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)

        cable2 = self.node.addChild("cable2")
        cable2.addObject("MechanicalObject", name='dof2', position=cableGeometry2)
        cable2.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry2))),
                        maxPositiveDisp='60',
                        maxDispVariation='1.5',
                        hasPullPoint='0',
                        minForce=0)
        cable2.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)

        cable3 = self.node.addChild("cable3")
        cable3.addObject("MechanicalObject", name='dof3', position=cableGeometry3)
        cable3.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry3))),
                        hasPullPoint='0',
                        maxPositiveDisp='60',
                        maxDispVariation='1.5',
                        minForce=0)
        cable3.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)
        
        cable4 = self.node.addChild("cable4")
        cable4.addObject("MechanicalObject", name='dof3', position=cableGeometry4)
        cable4.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry4))),
                        hasPullPoint='0',
                        maxPositiveDisp='130',
                        maxDispVariation='1.5',
                        minForce=0)
        cable4.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)
        
        cable5 = self.node.addChild("cable5")
        cable5.addObject("MechanicalObject", name='dof3', position=cableGeometry5)
        cable5.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry5))),
                        hasPullPoint='0',
                        maxPositiveDisp='130',
                        maxDispVariation='1.5',
                        minForce='0')
        cable5.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)
        
        cable6 = self.node.addChild("cable6")
        cable6.addObject("MechanicalObject", name='dof3', position=cableGeometry6)
        cable6.addObject('CableActuator', template='Vec3', name='CableActuator',
                        indices=list(range(len(cableGeometry6))),
                        hasPullPoint='0',
                        maxPositiveDisp='130',
                        maxDispVariation='1.5',
                        minForce='0')
        cable6.addObject('BarycentricMapping', name="Mapping", mapForces=False, mapMasses=False)  
        
    #Adding collision model to the your robot 
    def addCollisionModel(self, selfCollision=False):
        throatColli = self.node.addChild('CollisionModel')
        throatColli.addObject('MeshSTLLoader', name='loader', filename=path+'728_v1.stl')
        throatColli.addObject('MeshTopology', src='@loader')
        throatColli.addObject('MechanicalObject')
        #different collision model 
        #the models in different will have interaction
        throatColli.addObject('TriangleCollisionModel', contactStiffness=1, group=1,selfCollision=False)
        throatColli.addObject('LineCollisionModel' , contactStiffness=1, group=1,selfCollision=False)
        throatColli.addObject('PointCollisionModel', contactStiffness=1, group=1 ,selfCollision=False)
        throatColli.addObject('BarycentricMapping')
        
    #Adding visual model to the your robot 
    def addVisualModel(self, color):
        throatVisu = self.node.addChild('VisualModel')
        throatVisu.addObject('MeshSTLLoader', filename=path+'728_v1.stl')
        throatVisu.addObject('OglModel', color=color)
        throatVisu.addObject('BarycentricMapping')
    
    #Adding an effectors to your robot: It's the decision point 
    def addEffectors(self, target, position,rootNode):
        #Using drawing circle mode, just using this animation and activate the animate
        effectors = self.node.addChild('Effectors')
        effectors.addObject('MechanicalObject', position=position,showVectors=False,showVectorsScale=1)
        effectors.addObject('PositionEffector', indices=list(range(len(position))), effectorGoal=target_position)
        effectors.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
        effectors.addObject('Monitor', indices=0, showTrajectories=True,showPositions=False,TrajectoriesColor=[0,0,1,1],sizeFactor=10)
    def addTip(self):
        tip1=self.node.addChild('Tip1')
        tip1.addObject('MechanicalObject',name='tip1', position=[ 0, cablelength,130],)
        tip1.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
        tip2=self.node.addChild('Tip2')
        tip2.addObject('MechanicalObject',name='tip2', position=[ cablelength*cos(pi/6),-cablelength*cos(pi/3),130])
        tip2.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
        tip3=self.node.addChild('Tip3')
        tip3.addObject('MechanicalObject',name='tip3', position=[ -cablelength*cos(pi/6),-cablelength*cos(pi/3),130])
        tip3.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
    def addFixedbox(self):
        FixedBox=self.node.addChild("FixedBox")
        FixedBox.addObject('BoxROI', name='ROI1',box='5 5 0 -5 -5 -10', drawBoxes='true')
        FixedBox.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

        
##########################################
# Patient (Static or Dynamic obstacle )  #
##########################################
#Building model in different file, you need to do some rotation and translation
# rotation=[0,0,-90],translation=[0,10,150]

#patient-static without Dynamic obstacle
def patient_static(rootNode):
    patient=rootNode.addChild('patient')
    patient.addObject('MeshOBJLoader', name='loader', filename=path+'90tube-8in.obj', scale=1, rotation=[-90,90,0],translation=[0,0,180])
    patient.addObject('OglModel', src='@loader', color=[1,1,1,1])
    patient.addObject('MeshTopology', src='@loader', name='topo')
    patient.addObject('MechanicalObject')
    patient.addObject('TriangleCollisionModel', contactStiffness=1, group=2)
    patient.addObject('LineCollisionModel', contactStiffness=1, group=2)
    patient.addObject('PointCollisionModel', contactStiffness=1, group=2)
    return patient
def patient_real(rootNode):
    TranslationBase=[5850,-6190,6100]
    littleTongueBase=[TranslationBase[0]+3.5,TranslationBase[1]+13,TranslationBase[2]+150]
    patient_real=rootNode.addChild('patient_real')
    patient_real.addObject('MeshOBJLoader', name='loader', filename=path+'RealPatientWithoutGlottis.obj',
                           rotation=[0,180,0],translation=[littleTongueBase[0]+modelX,littleTongueBase[1]+modelY,littleTongueBase[2]+modelZ])
    
    # patient_real.addObject('MeshOBJLoader', name='loader', filename=path+'728final.obj',
    #                        rotation=[0,190,0],translation=[modelX,modelY,modelZ])
    patient_real.addObject('OglModel', src='@loader', color=[1,1,1,1])
    patient_real.addObject('MeshTopology', src='@loader', name='topo')
    patient_real.addObject('MechanicalObject')
    patient_real.addObject('TriangleCollisionModel', contactStiffness=1, group=2)
    patient_real.addObject('LineCollisionModel', contactStiffness=1, group=2)
    patient_real.addObject('PointCollisionModel', contactStiffness=1, group=2)
    return patient_real
def glottis(rootNode):
    TranslationBase=[5850,-6190,6100]
    littleTongueBase=[TranslationBase[0]+3.5,TranslationBase[1]+13,TranslationBase[2]+150]
    glottis=rootNode.addChild('glottis')
    glottis.addObject('MeshOBJLoader', name='loader', filename=path+'glottis.obj',
                      rotation=[0,180,0],translation=[littleTongueBase[0]+modelX,littleTongueBase[1]+modelY,littleTongueBase[2]+modelZ])
    # glottis.addObject('MeshOBJLoader', name='loader', filename=path+'glottis.obj',
    #                   rotation=[0,180,0],translation=[littleTongueBase[0],littleTongueBase[1],littleTongueBase[2]])
    glottis.addObject('OglModel', src='@loader', color=[1,1,1,1])
    return glottis
#patient with Dynamic obstacle 
#Using SurfacePressureActuator to simulate breathe state
def patient(rootNode):
    #add some solvers and load vtk file
    patient = rootNode.addChild('patient')
    patient.addObject('EulerImplicitSolver', name='odesolver', rayleighStiffness=0.1, rayleighMass=0.1)
    patient.addObject('SparseLDLSolver', name='preconditioner', template="CompressedRowSparseMatrixMat3x3d")
    patient.addObject('MeshVTKLoader', name='loader', filename=path+'tube_ball.vtk',rotation=[-90,90,0],translation=[0,0,210])
    patient.addObject('MeshTopology', src='@loader', name='container')
    patient.addObject('MechanicalObject', name='tetras', template='Vec3', showIndices=False, showIndicesScale=4e-5)
    patient.addObject('UniformMass', totalMass=0.5)
    patient.addObject('TetrahedronFEMForceField', template='Vec3', name='FEM', method='large', poissonRatio=0.3,  youngModulus=200)
    #Fix the robot in the environment
    patient.addObject('BoxROI', name='boxROI', box='-14 -10 148 14 12 178', drawBoxes=True)
    patient.addObject('RestShapeSpringsForceField', points='@boxROI.indices', stiffness=1e12)
    #Rigidization (necessary), the box need to include the whole mode
    patient.addObject('BoxROI', name='boxROISubTopo', box='-14 -45 148 14 14 225', drawBoxes=False)
    #This LinearSolverConstraintCorrection must be put in this position
    patient.addObject('LinearSolverConstraintCorrection', solverName='preconditioner')
    
    #Setting the goal which obstacle will get
    goal = rootNode.addChild('goal')
    goal.addObject('EulerImplicitSolver', firstOrder=True)
    goal.addObject('CGLinearSolver', iterations='100', tolerance='1e-05', threshold='1e-05')
    goal.addObject('MechanicalObject', name='goalMO', position='0 -10 205',showObject=True, showObjectScale=2)
    goal.addObject('SphereCollisionModel', radius='0.25', group='3')
    goal.addObject('UncoupledConstraintCorrection')
    
    #Animate goal.position and it will move
    def myAnimation(target, factor):  # here you define your anima # factor = (currentTime-startTime) / duration
        z = 205-5*factor
        target.value = [[0, -10,z]]
    animate(myAnimation, {'target': goal.goalMO.position}, duration=0.5,mode="pingpong")
    
    #Create a cavity for the SurfacePressureActuator methods
    cavity1 = patient.addChild('cavity1')
    cavity1.addObject('MeshSTLLoader', name='loader', filename=path+'cavity1.stl',rotation=[-90,90,0],translation=[0,0,210])
    cavity1.addObject('MeshTopology', src='@loader', name='topo')
    cavity1.addObject('MechanicalObject', name='cavity')
    cavity1.addObject('SurfacePressureActuator', template='Vec3', triangles='@topo.triangles', minPressure=0, maxVolumeGrowthVariation=500)
    cavity1.addObject('BarycentricMapping', name='mapping',  mapForces=False, mapMasses=False)
    
    #Fix part of the cavity
    patient.addObject('BoxROI', name='boxROIball1', drawBoxes=True,box=[12, -22, 225, -12, 2, 210])
    patient.addObject('RestShapeSpringsForceField', points='@boxROIball1.indices', stiffness=1e12, angularStiffness=1e12)

    #Add an effector on patient 
    effector = patient.addChild('effector')
    effector.addObject('MechanicalObject', name="effectorPoint", position="0 -10 205",showVectorsScale=1)
    effector.addObject('PositionEffector', template='Vec3', indices=0, effectorGoal="@../../goal/goalMO.position")
    # effector.addObject('PositionEffector', template='Vec3', indices=0, effectorGoal=target)
    effector.addObject('BarycentricMapping', mapForces=False, mapMasses=False)
    
    #Add a visual model
    patientVisu = patient.addChild('visu')
    patientVisu.addObject('MeshObjLoader', filename=path+"tube_ball.stl", name="loader",rotation=[-90,90,0],translation=[0,0,210])
    patientVisu.addObject('OglModel', src="@loader", template='Vec3', color=[0.7, 0.7, 0.7, 0.2])
    patientVisu.addObject('BarycentricMapping')
    
    #Add collision model to patient
    collisionPatient = patient.addChild('collisionPatient')
    collisionPatient.addObject('MeshSTLLoader', name='loader', filename=path+'tube_ball.stl', rotation=[-90,90,0],translation=[0,0,210])
    collisionPatient.addObject('MeshTopology', src='@loader', name='topo')
    collisionPatient.addObject('MechanicalObject', name='collisMech')
    collisionPatient.addObject('TriangleCollisionModel', selfCollision=False, contactStiffness=1, group=2)
    collisionPatient.addObject('LineCollisionModel', selfCollision=False, contactStiffness=1, group=2)
    collisionPatient.addObject('PointCollisionModel', selfCollision=False, contactStiffness=1, group=2)
    collisionPatient.addObject('BarycentricMapping')

#Using _regidify function to set your robot into deformable and rigid parts     

def createScene(rootNode):
    #RequiredPlugins
    rootNode.addObject('RequiredPlugin', pluginName=['SoftRobots','SofaValidation','Sofa.GL.Component.Rendering2D','CImgPlugin','Sofa.Component.SolidMechanics.Spring','SofaMeshCollision','ArticulatedSystemPlugin','Sofa.Component.Mapping.MappedMatrix','Sofa.Component.Topology.Container.Dynamic','SofaSparseSolver','SofaPreconditioner','SofaPython3','SofaConstraint',
                                                     'SofaImplicitOdeSolver','SofaLoader','SofaSimpleFem','SofaBoundaryCondition','SofaEngine',
                                                     'SofaOpenglVisual','Sofa.Component.IO.Mesh','Sofa.Component.LinearSolver.Direct','Sofa.Component.LinearSolver.Iterative',
                                                     'Sofa.Component.Mass','Sofa.Component.ODESolver.Backward','Sofa.Component.SceneUtility',
                                                     'Sofa.Component.StateContainer','Sofa.Component.Topology.Container.Constant','Sofa.Component.Visual',
                                                     'Sofa.GL.Component.Rendering3D','Sofa.Component.SolidMechanics.FEM.Elastic','Sofa.Component.AnimationLoop',
                                                     'Sofa.Component.Collision.Detection.Algorithm','Sofa.Component.Collision.Detection.Intersection',
                                                     'Sofa.Component.Collision.Geometry','Sofa.Component.Collision.Response.Contact','Sofa.Component.Constraint.Lagrangian.Correction',
                                                     'Sofa.Component.Constraint.Projective','Sofa.Component.Mapping.Linear','Sofa.Component.Engine.Select',
                                                     'Sofa.Component.Setting'])
    rootNode.addObject('RequiredPlugin', name='SoftRobots.Inverse')
    
    rootNode.addObject('DefaultVisualManagerLoop')
    #If you want show anything like vtk model, visual model 
    # rootNode.addObject('VisualStyle', displayFlags='showBehavior')
    #A clear visual help you observe
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideBoundingCollisionModels hideForceFields showInteractionForceFields hideWireframe')
    
    #Set gravity of your scene
    rootNode.gravity = [0., 0., 9.810]
    # rootNode.gravity = [0., 0., 0.]
    
    #Add  AnimationManager  for using myAnimation function
    rootNode.addObject(AnimationManager(rootNode))
    rootNode.addObject('FreeMotionAnimationLoop')
    
    #Set the background, also you can add a image file 
    # rootNode.addObject('BackgroundSetting',color=[0,0,0,1])
    
    #QPInverseProblemSolver for inverse mode
    rootNode.addObject('QPInverseProblemSolver',epsilon=1e-1)
    # rootNode.addObject('SofaSparseSolver')

    #Add a new child node for your throat 
    simulation = rootNode.addChild('Simulation')
    #Necessary things
    simulation.addObject('EulerImplicitSolver', name='odesolver', firstOrder=False, rayleighMass=0.1, rayleighStiffness=0.1)
    simulation.addObject('ShewchukPCGLinearSolver', name='linearSolver', iterations=500, tolerance=1.0e-18, preconditioners='precond')
    simulation.addObject('SparseLDLSolver',template="CompressedRowSparseMatrixMat3x3d", name='precond')
    simulation.addObject('GenericConstraintCorrection', solverName='precond')

    #Define your robot, add visual model and collision model
    throat = Throat(simulation)
    throat.addVisualModel(color=[1., 1., 1., 0.8])
    throat.addCollisionModel()
    
    
    
    #Make some points bigger, this step is not necessary
    #It only help to observe points on each parts
    # setData(simulation.RigidifiedStructure.RigidParts.dofs, showObject=True, showObjectScale=1, drawMode=2)
    # setData(simulation.RigidifiedStructure.RigidParts.RigidifiedParticules.dofs, showObject=True, showObjectScale=0.8,drawMode=1, showColor=[1., 1., 1., 1.])
    

    ##########################################
    # Select your patient                    #
    ##########################################
    # patient(rootNode)
    # patient_static(rootNode)
    # patient_real(rootNode)
    # glottis(rootNode)
    #Open interaction for collision model
    rootNode.addObject('DefaultPipeline')
    rootNode.addObject('BruteForceBroadPhase')
    rootNode.addObject('BVHNarrowPhase')
    rootNode.addObject('DefaultContactManager', response="FrictionContact", responseParams="mu=0")
    rootNode.addObject('LocalMinDistance', alarmDistance=2, contactDistance=0.5)
    
    #Set a target for your robot and activate your effector on robot
    # effectorTarget2(rootNode)
    target = effectorTarget(rootNode)
    throat.addEffectors(target=target.dofs.getData('position').getLinkPath(), position=robot_decisionpoint,rootNode=rootNode)
    
    ##########################################
    #Clip the whole model                    #
    ##########################################
    # rootNode.addObject("ClipPlane",name="Clip" ,normal="1 0 0")
    
    ##########################################
    #Camera                                  #
    ##########################################
    myValueFile=[]
    myForcesFile=[]
    def myAnimation(tip1,tip2,tip3,camera_position1,cable1,cable2,cable3,cable4,cable5,cable6,print_flag,factor,effector,target,rootNode,hwnd_title):
        #Calculate camera position and orientation
        x1=tip1.tip1.position.value[0][0]
        y1=tip1.tip1.position.value[0][1]
        z1=tip1.tip1.position.value[0][2]
        x2=tip2.tip2.position.value[0][0]
        y2=tip2.tip2.position.value[0][1]
        z2=tip2.tip2.position.value[0][2]
        x3=tip3.tip3.position.value[0][0]
        y3=tip3.tip3.position.value[0][1]
        z3=tip3.tip3.position.value[0][2]
        p1=np.array([x1,y1,z1])
        p2=np.array([x2,y2,z2])
        p3=np.array([x3,y3,z3])
        A=p1-p2
        B=p3-p2
        A_B=np.cross(A,B)
        x=A_B[0]
        y=A_B[1]
        z=A_B[2]
        angleX=acos(sqrt((x**2 + z**2)/(x**2+y**2+z**2)))
        angleY=atan(x/z)
        angleZ=0
        q = Quat.createFromEuler([angleX, angleY, angleZ], 'sxyz')
        ret1=R.from_quat(q)
        matrix=ret1.as_matrix()
        rot=np.mat([[-1,0,0],[0,1,0],[0,0,-1]])
        matrix_rot=matrix*rot
        ret3=R.from_matrix(matrix_rot)
        q_final=ret3.as_quat()
        # print(q_final)
        camera_position1.position=[(x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3]
        camera_position1.orientation=[q_final[0],q_final[1],q_final[2],q_final[3]]
        #Save Data
        # print(effector.PositionEffector.effectorGoal.value)
        node=rootNode
        simNode=node.getChild('Simulation')
        throatNode=simNode.getChild('Throat')
        effNode=throatNode.getChild('Effectors')
        print(rootNode.time.value)
        # print(target.position.value[0][0])
        if target.position.value[0][0]+5>effector.MechanicalObject.position.value[0][0]>target.position.value[0][0]-5 and target.position.value[0][1]+5>effector.MechanicalObject.position.value[0][1]>target.position.value[0][1]-5 and target.position.value[0][2]+5>effector.MechanicalObject.position.value[0][2]>target.position.value[0][2]-5:      
            newx=effNode.MechanicalObject.position.value[0][0]
            newy=effNode.MechanicalObject.position.value[0][1]
            newz=effNode.MechanicalObject.position.value[0][2]
            newGoal=[newx,newy,newz]
            print('Mission Completed! The newGoal is'+str(newGoal))
            
            effNode.PositionEffector.effectorGoal.value=[newGoal]
            # effNode.PositionEffector.effectorGoal.value=[target.position.value[0]]
            myValueFile.append([cable1.CableActuator.displacement.value,cable2.CableActuator.displacement.value,cable3.CableActuator.displacement.value,
                                cable4.CableActuator.displacement.value,cable5.CableActuator.displacement.value,cable6.CableActuator.displacement.value,
                                factor,
                                effector.MechanicalObject.position.value[0][0],effector.MechanicalObject.position.value[0][1],effector.MechanicalObject.position.value[0][2]])
            myForcesFile.append([cable1.CableActuator.force.value,cable2.CableActuator.force.value,cable3.CableActuator.force.value,
                                cable4.CableActuator.force.value,cable5.CableActuator.force.value,cable6.CableActuator.force.value,
                                factor])
        else:
            effNode.PositionEffector.effectorGoal.value=[target.position.value[0]]
            myValueFile.append([cable1.CableActuator.displacement.value,cable2.CableActuator.displacement.value,cable3.CableActuator.displacement.value,
                                cable4.CableActuator.displacement.value,cable5.CableActuator.displacement.value,cable6.CableActuator.displacement.value,
                                factor,
                                effector.MechanicalObject.position.value[0][0],effector.MechanicalObject.position.value[0][1],effector.MechanicalObject.position.value[0][2]])
            print('Get Value!'+str(factor))
            myForcesFile.append([cable1.CableActuator.force.value,cable2.CableActuator.force.value,cable3.CableActuator.force.value,
                                cable4.CableActuator.force.value,cable5.CableActuator.force.value,cable6.CableActuator.force.value,
                                factor])

            
        def get_all_hwnd(hwnd,mouse):
            if win32gui.IsWindow(hwnd) and win32gui.IsWindowEnabled(hwnd) and win32gui.IsWindowVisible(hwnd):
                hwnd_title.update({hwnd:win32gui.GetWindowText(hwnd)})

        # win32gui.EnumWindows(get_all_hwnd, 0)
        # find the window name
        # for h,t in hwnd_title.items():
        #     if t != "":
        #         print(h, t)
        # # ScreenShot
        # hwnd = win32gui.FindWindow(None, 'SOFA v22.06.99 - C:/Users/TeoRen/Desktop/sim2real_ThroatRobotwithPatientlocal/throat_v60_InverseMode.py')
        # screen = QApplication.primaryScreen()
        # img = screen.grabWindow(hwnd).toImage()
        # img.save(Img_path+str(factor)+"_screenshot.jpg")    
        
    def myOnDone1(tip1,tip2,tip3,camera_position1,cable1,cable2,cable3,cable4,cable5,cable6,print_flag,factor,effector,target,rootNode,hwnd_title):
        # pass
        # ser.close()
        
        if print_flag==0:
            # print(myValueFile)
            file = open(dirPath+'data/myValue.csv','a')
            for i in range(len(myValueFile)):
                s = str(myValueFile[i]).replace('[','').replace(']','')+'\n'
                file.write(s)
            file.close()
            
            file = open(dirPath+'data/myForcesValue.csv','a')
            for i in range(len(myForcesFile)):
                s = str(myForcesFile[i]).replace('[','').replace(']','')+'\n'
                file.write(s)
            file.close()
            print_flag+=1
        else:
            print("Aleady save the data!")
    
    # rootNode.addObject(Controller(rootNode))
    # rootNode.addObject("OglViewport", name="OglViewport",screenPosition="300 0",screenSize="400 400",zFar=140,zNear=1)
    rootNode.addObject("RecordedCamera",zFar=140,zNear=1)
    
    # animate(myAnimation, {'target': effectors.MechanicalObject,'camera_position1':rootNode.OglViewport,'camera_position2':rootNode.RecordedCamera,'cable4':self.node.cable4,'cable5':self.node.cable5,'cable6':self.node.cable6},duration=5)
    animate(myAnimation, {'tip1':rootNode.Simulation.Throat.Tip1,'tip2':rootNode.Simulation.Throat.Tip2,'tip3':rootNode.Simulation.Throat.Tip3,
                          'camera_position1':rootNode.RecordedCamera,
                          'cable1':rootNode.Simulation.Throat.cable1,'cable2':rootNode.Simulation.Throat.cable2,'cable3':rootNode.Simulation.Throat.cable3,
                          'cable4':rootNode.Simulation.Throat.cable4,'cable5':rootNode.Simulation.Throat.cable5,'cable6':rootNode.Simulation.Throat.cable6,
                          'print_flag':print_flag,
                          'effector':rootNode.Simulation.Throat.Effectors,
                          'target':rootNode.Target.dofs,'rootNode':rootNode,'hwnd_title':hwnd_title},duration=2)