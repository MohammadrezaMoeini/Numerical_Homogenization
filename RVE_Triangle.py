# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 11:29:43 2022
@author: momoe

This scrit creates a traingle RVE.
"""


# *****************************************************************************
' -------------------------------  Geometry -----------------------------------' 
# *****************************************************************************
'Defult Geometry'
import numpy as np

' Abaqus modules '
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
from part import *
from material import *
from section import *
from optimization import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import time
'''please dont delete the import lines in each functions, I tried to write just one function
instead of them and just execute one function but it is really better to write all the lines'''
# ---------------------------------
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import optimization
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import os


def Cube(X,Y,Z):
    """Designing the cube
    Parameters:
        X,Y,Z: half of the each side of the cube
    """ 
      
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2*int(Z))
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
        
    s.Line(point1=(-X, -Y), point2=(X, -Y))
    s.HorizontalConstraint(entity=g[2], addUndoState=False)
    
    s.Line(point1=(X, -Y), point2=(X, Y))
    s.VerticalConstraint(entity=g[3], addUndoState=False)    
    s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    
    s.Line(point1=(X, Y), point2=(-X, Y))    
    s.HorizontalConstraint(entity=g[4], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    
    s.Line(point1=(-X, Y), point2=(-X, -Y))
    s.VerticalConstraint(entity=g[5], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
    
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseSolidExtrude(sketch=s, depth=2.*Z)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

    
def compute_points(a, t):
    """
    Returns the required points to cut the internal and external traingles
        Parameters:
        a: the length size of the traingle.
        t: wall thickness in the lattice (Here, in the RVE, we need half of it)
    
        Note: It will be over a cube with dimension of (2X*2Y*2Z)
    """
    #------------------------- Right Side
    # Right_Internal traingle 
    p1 = np.array([+t/2.0, +a/2.0]) #ok
    p2 = np.array([+t/2. + a*np.sqrt(3.)/2., 0.0]) #ok
    p3 = np.array([+t/2.0, -a/2.0]) #ok
    # Top_Right traingle 
    p4 = np.array([+t/2. + t/2., t*np.sqrt(3.)/2. + a/2.]) #ok
    p5 = np.array([+t/2. + t/2. + a*np.sqrt(3.)/2., t*np.sqrt(3.)/2. + a/2.])
    p6 = np.array([+t/2. + t/2. + a*np.sqrt(3.)/2., t*np.sqrt(3.)/2.])
    # Bottom_Right traingle
    p7 = np.array([p4[0], -p4[1]])
    p8 = np.array([p6[0], -p6[1]])
    p9 = np.array([p5[0], -p5[1]])
    #------------------------- Left Side
    # Left_Internal traingle
    q1 = np.array([-p1[0], p1[1]])
    q2 = np.array([-p2[0], p2[1]])
    q3 = np.array([-p3[0], p3[1]])
    # Top_Left traingle
    q4 = np.array([-p4[0], p4[1]])
    q5 = np.array([-p5[0], p5[1]])
    q6 = np.array([-p6[0], p6[1]])
    # Bottom_Left traingle
    q7 = np.array([-p7[0], p7[1]])
    q8 = np.array([-p8[0], p8[1]])
    q9 = np.array([-p9[0], p9[1]])
    # Corners
    pp5 = np.array([p5[0]+t/2., p5[1]])
    pp9 = np.array([p9[0]+t/2., p9[1]])
    qq5 = np.array([q5[0]-t/2., q5[1]])
    qq9 = np.array([q9[0]-t/2., q9[1]])
    
    X = p5[0] + t/2.
    Z = p5[1]

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, q1, q2, q3, q4, q5, q6, q7, q8, q9, X, Z, pp5, pp9, qq5, qq9


def Cut_traingles(X, Y, Z,
               p1, p2, p3, p4, p5, p6, p7, p8, p9,
               q1, q2, q3, q4, q5, q6, q7, q8, q9):
    """Cutting the internal and external triangles on the top face of a cube
    Parameters:
        pi = the required points to cut the triangles at the right side. 
        qi = the required points to cut the triangles at the left side. 
    """

    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    t = p.MakeSketchTransform(sketchPlane=f1[0], sketchUpEdge=e[1], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=t)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    

    #------------------------- Right Side
    s1.Line(point1=(p1[0], p1[1]), point2=(p2[0], p2[1]))   
    s1.Line(point1=(p2[0], p2[1]), point2=(p3[0], p3[1]))
    s1.Line(point1=(p3[0], p3[1]), point2=(p1[0], p1[1]))
    
    s1.Line(point1=(p4[0], p4[1]), point2=(p5[0], p5[1]))
    s1.Line(point1=(p5[0], p5[1]), point2=(p6[0], p6[1]))
    s1.Line(point1=(p6[0], p6[1]), point2=(p4[0], p4[1]))
    
    s1.Line(point1=(p7[0], p7[1]), point2=(p8[0], p8[1]))
    s1.Line(point1=(p8[0], p8[1]), point2=(p9[0], p9[1]))
    s1.Line(point1=(p7[0], p9[1]), point2=(p9[0], p9[1]))
    
    #------------------------- Lef Side
    s1.Line(point1=(q1[0], q1[1]), point2=(q2[0], q2[1]))   
    s1.Line(point1=(q2[0], q2[1]), point2=(q3[0], q3[1]))
    s1.Line(point1=(q3[0], q3[1]), point2=(q1[0], q1[1]))
    
    s1.Line(point1=(q4[0], q4[1]), point2=(q5[0], q5[1]))
    s1.Line(point1=(q5[0], q5[1]), point2=(q6[0], q6[1]))
    s1.Line(point1=(q6[0], q6[1]), point2=(q4[0], q4[1]))
    
    s1.Line(point1=(q7[0], q7[1]), point2=(q8[0], q8[1]))
    s1.Line(point1=(q8[0], q8[1]), point2=(q9[0], q9[1]))
    s1.Line(point1=(q7[0], q7[1]), point2=(q9[0], q9[1]))

    
    p = mdb.models['Model-1'].parts['Part-1']
    f, e1 = p.faces, p.edges
    p.CutExtrude(sketchPlane=f[0], sketchUpEdge=e1[1], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=2.*Y, 
        flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']


def Fill_traingles(X, Y, Z,
                   p1, p2, p3, p4, p5, p6, p7, p8, p9,
                   q1, q2, q3, q4, q5, q6, q7, q8, q9):
    """Filling triangles on cube (to assing the elastic air)
    Parameters:
        pi = the required points to cut the triangles at the right side. 
        qi = the required points to cut the triangles at the left side. 
    """
    
    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    tp = p.MakeSketchTransform(sketchPlane=f1[19], sketchUpEdge=e[61], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=tp)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    #------------------------- Right Side
    s1.Line(point1=(p1[0], p1[1]), point2=(p2[0], p2[1]))   
    s1.Line(point1=(p2[0], p2[1]), point2=(p3[0], p3[1]))
    s1.Line(point1=(p3[0], p3[1]), point2=(p1[0], p1[1]))
    
    s1.Line(point1=(p4[0], p4[1]), point2=(p5[0], p5[1]))
    s1.Line(point1=(p5[0], p5[1]), point2=(p6[0], p6[1]))
    s1.Line(point1=(p6[0], p6[1]), point2=(p4[0], p4[1]))
    
    s1.Line(point1=(p7[0], p7[1]), point2=(p8[0], p8[1]))
    s1.Line(point1=(p8[0], p8[1]), point2=(p9[0], p9[1]))
    s1.Line(point1=(p7[0], p9[1]), point2=(p9[0], p9[1]))
    
    #------------------------- Lef Side
    s1.Line(point1=(q1[0], q1[1]), point2=(q2[0], q2[1]))   
    s1.Line(point1=(q2[0], q2[1]), point2=(q3[0], q3[1]))
    s1.Line(point1=(q3[0], q3[1]), point2=(q1[0], q1[1]))
    
    s1.Line(point1=(q4[0], q4[1]), point2=(q5[0], q5[1]))
    s1.Line(point1=(q5[0], q5[1]), point2=(q6[0], q6[1]))
    s1.Line(point1=(q6[0], q6[1]), point2=(q4[0], q4[1]))
    
    s1.Line(point1=(q7[0], q7[1]), point2=(q8[0], q8[1]))
    s1.Line(point1=(q8[0], q8[1]), point2=(q9[0], q9[1]))
    s1.Line(point1=(q7[0], q7[1]), point2=(q9[0], q9[1]))

    p = mdb.models['Model-1'].parts['Part-1']
    f1, e1 = p.faces, p.edges
    p.SolidExtrude(sketchPlane=f1[19], sketchUpEdge=e1[61], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=2.*Y, 
        flipExtrudeDirection=ON, keepInternalBoundaries=ON)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']


# *****************************************************************************
' -------------------------------  Material -----------------------------------' 
# *****************************************************************************
def Material_Section_Solid(E, v, Material_name, Sections_name):

    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    mdb.models['Model-1'].Material(name=Material_name)
    mdb.models['Model-1'].materials[Material_name].Elastic(table=((E, v), 
        ))
    mdb.models['Model-1'].HomogeneousSolidSection(name= Sections_name, 
        material=Material_name, thickness=None)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(cells=cells, name='Set-1')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName= Sections_name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


def Material_Section_Void(E, v, Material_name, Sections_name):

    mdb.models['Model-1'].Material(name=Material_name)
    mdb.models['Model-1'].materials[Material_name].Elastic(table=((E, v), ))
    mdb.models['Model-1'].HomogeneousSolidSection(name=Sections_name, material=Material_name, 
        thickness=None)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#7e ]', ), )
    region = p.Set(cells=cells, name='Set-2')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName=Sections_name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

def Material_color():

    session.viewports['Viewport: 1'].enableMultipleColors()
    session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
    cmap = session.viewports['Viewport: 1'].colorMappings['Material']
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors()
    session.viewports['Viewport: 1'].enableMultipleColors()
    session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
    cmap = session.viewports['Viewport: 1'].colorMappings['Material']
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors()

# *****************************************************************************
' -----------------------------  Assembly -------------------------------------' 
# ***************************************************************************** 
def Assembly():
    """ABAQUS_function: Assembly"""

    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name=
    'Part-1', part=
    mdb.models['Model-1'].parts['Part-1'])

def Move_Ref_to(Xp, Yp, Zp):
    """ABAQUS_function:
        Move reference point (0,0,0) to (Xp,Yp,Zp)"""
        
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1', ), vector=
    (Xp, Yp, Zp))



# *****************************************************************************
' -------------------------------  Step --------------------------------------' 
# *****************************************************************************    
    
def Step():
    """ABAQUS_function: Step:
           Return: Stress, Strain, Elastic energy """
   
           
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'E', 'U', 'IVOL', 'ELSE'))


# *****************************************************************************
' -------------------------------  Mesh  --------------------------------------' 
# ***************************************************************************** 
def Mesh(element_size,min_size_fac, dF):
    """
    Element type: Hex dominated element
    Algorithm: advanced front
    Order: Quadratic element 
    """
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedRegions = c1.getSequenceFromMask(mask=('[#7f ]', ), )
    a.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED, 
        technique=SWEEP)
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.seedPartInstance(regions=partInstances, 
                       size=element_size, 
                       deviationFactor=dF, 
        minSizeFactor=min_size_fac)
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.generateMesh(regions=partInstances)
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    cells1 = c1.getSequenceFromMask(mask=('[#7f ]', ), )
    pickedRegions =(cells1, )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))

