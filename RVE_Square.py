# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:24:13 2022
@author: momoe

This script generates RVE for a square lattice structure.
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
    


def Cut_square(X, Y, Z, p1, p2, p3, p4):
    """Cutting the internal square on the top face of a cube
    Parameters:
        pi = the required points to cut the internal square. 
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
    

    
    s1.Line(point1=(p1[0], p1[1]), point2=(p2[0], p2[1]))   
    s1.Line(point1=(p2[0], p2[1]), point2=(p3[0], p3[1]))   
    s1.Line(point1=(p3[0], p3[1]), point2=(p4[0], p4[1]))       
    s1.Line(point1=(p4[0], p4[1]), point2=(p1[0], p1[1]))
    
    p = mdb.models['Model-1'].parts['Part-1']
    f, e1 = p.faces, p.edges
    p.CutExtrude(sketchPlane=f[0], sketchUpEdge=e1[1], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=2.*Y, 
        flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    
    
def compute_points(X, Y, Z, t):
    """
    Returns the required points to cut the internal square.
        Parameters:
        t: wall thickness in the lattice (Here, in the RVE, we need half of it
        as tt)
    """
    tt = t/2.0
    p1 = np.array([X - tt, Y - tt])
    p2 = np.array([-X + tt, Y - tt])
    p3 = np.array([-X + tt, -Y + tt])
    p4 = np.array([X - tt, -Y + tt])
    
    return p1, p2, p3, p4
    


def Fill_square(X, Y, Z, t, p1, p2, p3, p4):
    """Filling hexagonal on cube (to assing the elastic air)
    Parameters:
        pi: 6 points of the hexgonal.
    """
    
    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    tp = p.MakeSketchTransform(sketchPlane=f1[4], sketchUpEdge=e[13], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=tp)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    s1.Line(point1=(p1[0], p1[1]), point2=(p2[0], p2[1]))   
    s1.Line(point1=(p2[0], p2[1]), point2=(p3[0], p3[1]))   
    s1.Line(point1=(p3[0], p3[1]), point2=(p4[0], p4[1]))       
    s1.Line(point1=(p4[0], p4[1]), point2=(p1[0], p1[1]))

    p = mdb.models['Model-1'].parts['Part-1']
    f1, e1 = p.faces, p.edges
    p.SolidExtrude(sketchPlane=f1[4], sketchUpEdge=e1[13], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=2.*Y, 
        flipExtrudeDirection=ON, keepInternalBoundaries=ON)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    

# *****************************************************************************
' -------------------------------  Material -------------------------------------' 
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
    cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
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
' -----------------------------  Assembly ------------------------------------' 
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
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.seedPartInstance(regions=partInstances, 
                       size = element_size,
                       deviationFactor = min_size_fac, 
                       minSizeFactor = dF)
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.generateMesh(regions=partInstances)


def make_it_Quadratic():
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    cells1 = c1.getSequenceFromMask(mask=('[#3 ]', ), )
    pickedRegions =(cells1, )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))


