# ***************************************************************************************************
' -------------------------------  Geometry -------------------------------------' 
# ***************************************************************************************************
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
'''please do not delete the import lines in each functions, I tried to write just one function
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
    
    
def Cut_hex(h, Y, Z, p1, p2, p3, p4, p5, p6):
    """Cutting hexagonal on the top face of a cube
    Parameters:
        pi: 6 points of the hexgonal.
        h: Depth of the cut
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
    s1.Line(point1=(p4[0], p4[1]), point2=(p5[0], p5[1]))
    s1.Line(point1=(p5[0], p5[1]), point2=(p6[0], p6[1]))   
    s1.Line(point1=(p6[0], p6[1]), point2=(p1[0], p1[1]))
    
    p = mdb.models['Model-1'].parts['Part-1']
    f, e1 = p.faces, p.edges
    p.CutExtrude(sketchPlane=f[0], sketchUpEdge=e1[1], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=2.*Y, 
        flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
   
    
def Cut_hex_external(h, Y, Z, q1, q2, q3, q4):
    """ Cutting external part on the top face of cuve"""

    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    t = p.MakeSketchTransform(sketchPlane=f1[6], sketchUpEdge=e[19], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=t)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    s1.Line(point1=(q1[0], q1[1]), point2=(q2[0], q2[1]))    
    s1.Line(point1=(q2[0], q2[1]), point2=(q3[0], q3[1]))    
    s1.Line(point1=(q3[0], q3[1]), point2=(q4[0], q4[1]))        
    s1.Line(point1=(q4[0], q4[1]), point2=(q1[0], q1[1]))
        
    s1.Line(point1=(-q1[0], q1[1]), point2=(-q2[0], q2[1]))    
    s1.Line(point1=(-q2[0], q2[1]), point2=(-q3[0], q3[1]))    
    s1.Line(point1=(-q3[0], q3[1]), point2=(-q4[0], q4[1]))        
    s1.Line(point1=(-q4[0], q4[1]), point2=(-q1[0], q1[1]))
    
    s1.Line(point1=(-q1[0], -q1[1]), point2=(-q2[0], -q2[1]))    
    s1.Line(point1=(-q2[0], -q2[1]), point2=(-q3[0], -q3[1]))    
    s1.Line(point1=(-q3[0], -q3[1]), point2=(-q4[0], -q4[1]))        
    s1.Line(point1=(-q4[0], -q4[1]), point2=(-q1[0], -q1[1]))
    
    s1.Line(point1=(q1[0], -q1[1]), point2=(q2[0], -q2[1]))    
    s1.Line(point1=(q2[0], -q2[1]), point2=(q3[0], -q3[1]))    
    s1.Line(point1=(q3[0], -q3[1]), point2=(q4[0], -q4[1]))        
    s1.Line(point1=(q4[0], -q4[1]), point2=(q1[0], -q1[1]))

    p = mdb.models['Model-1'].parts['Part-1']
    f, e1 = p.faces, p.edges
    p.CutExtrude(sketchPlane=f[6], sketchUpEdge=e1[19], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=h, 
        flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']


def Fill_hex(h, Y, Z, t, p1, p2, p3, p4, p5, p6):
    """Cutting hexagonal on the top face of a cube
    Parameters:
        pi: 6 points of the hexgonal.
        h: Depth of the cut
    """
    
    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    tp = p.MakeSketchTransform(sketchPlane=f1[10], sketchUpEdge=e[31], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, t-Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=tp)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    s1.Line(point1=(p1[0], p1[1]), point2=(p2[0], p2[1]))   
    s1.Line(point1=(p2[0], p2[1]), point2=(p3[0], p3[1]))   
    s1.Line(point1=(p3[0], p3[1]), point2=(p4[0], p4[1]))       
    s1.Line(point1=(p4[0], p4[1]), point2=(p5[0], p5[1]))
    s1.Line(point1=(p5[0], p5[1]), point2=(p6[0], p6[1]))   
    s1.Line(point1=(p6[0], p6[1]), point2=(p1[0], p1[1]))
    
    p = mdb.models['Model-1'].parts['Part-1']
    f1, e1 = p.faces, p.edges
    p.SolidExtrude(sketchPlane=f1[10], sketchUpEdge=e1[31], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=h, 
        flipExtrudeDirection=ON, keepInternalBoundaries=ON)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    

    
def Fill_hex_external(h, Y, Z, t, q1, q2, q3, q4):
    """Cutting hexagonal on the top face of a cube
    Parameters:
        pi: 6 points of the hexgonal.
        h: Depth of the cut
    """

    p = mdb.models['Model-1'].parts['Part-1']
    f1, e = p.faces, p.edges
    
    tp = p.MakeSketchTransform(sketchPlane=f1[12], sketchUpEdge=e[43], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, t-Y, Z))
    
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=2*int(Z), gridSpacing=0.2, transform=tp)
    
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p = mdb.models['Model-1'].parts['Part-1']
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    s1.Line(point1=(q1[0], q1[1]), point2=(q2[0], q2[1]))    
    s1.Line(point1=(q2[0], q2[1]), point2=(q3[0], q3[1]))    
    s1.Line(point1=(q3[0], q3[1]), point2=(q4[0], q4[1]))        
    s1.Line(point1=(q4[0], q4[1]), point2=(q1[0], q1[1]))
        
    s1.Line(point1=(-q1[0], q1[1]), point2=(-q2[0], q2[1]))    
    s1.Line(point1=(-q2[0], q2[1]), point2=(-q3[0], q3[1]))    
    s1.Line(point1=(-q3[0], q3[1]), point2=(-q4[0], q4[1]))        
    s1.Line(point1=(-q4[0], q4[1]), point2=(-q1[0], q1[1]))
    
    s1.Line(point1=(-q1[0], -q1[1]), point2=(-q2[0], -q2[1]))    
    s1.Line(point1=(-q2[0], -q2[1]), point2=(-q3[0], -q3[1]))    
    s1.Line(point1=(-q3[0], -q3[1]), point2=(-q4[0], -q4[1]))        
    s1.Line(point1=(-q4[0], -q4[1]), point2=(-q1[0], -q1[1]))
    
    s1.Line(point1=(q1[0], -q1[1]), point2=(q2[0], -q2[1]))    
    s1.Line(point1=(q2[0], -q2[1]), point2=(q3[0], -q3[1]))    
    s1.Line(point1=(q3[0], -q3[1]), point2=(q4[0], -q4[1]))        
    s1.Line(point1=(q4[0], -q4[1]), point2=(q1[0], -q1[1]))
    
    p = mdb.models['Model-1'].parts['Part-1']
    f1, e1 = p.faces, p.edges
    p.SolidExtrude(sketchPlane=f1[12], sketchUpEdge=e1[43], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s1, depth=h, 
        flipExtrudeDirection=ON, keepInternalBoundaries=ON)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
        
    
       
    
# .................
# # Geometrical Parameters ---> calculator
# ...............
# =============================================================================
# Dimentions of the cube WRT c,t,h
# =============================================================================
def f_cth_to_XZY(c,t,h):
    '''independents variables c,t,h to---> X,Y,Z'''
    X=3.*c + np.sqrt(3)*t
    Z=np.sqrt(3)*c + t
    Y=0.5*(h)
    
    return X,Z,Y
# =============================================================================
# Points of the hexagonale 
# =============================================================================
def hexagonale_points(c,t,X,Z):
    """Calculator: points of hexagonal WRT X & Z"""
    
    p1=[c, Z-t]
    p2=[2.*c, 0.]
    
    p3=[p1[0], -p1[1]]    
    p4=[-p1[0], -p1[1]]
    
    p5=[-p2[0], p2[1]]    
    p6=[-p1[0], p1[1]]
    
    return p1, p2, p3, p4, p5, p6
# =============================================================================
# Points of the hexagonale (EXTERNAL)
# =============================================================================
def hexagonale_points_External(c,t,X,Z):
    """Calculator: points of external points WRT X & Z"""
    
    q1=[c+(np.sqrt(3.))*t, Z]
    q2=[X-c, t]
    q3=[X, t]
    q4=[X,Z]
        
    return q1,q2,q3,q4

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
    cells = c.getSequenceFromMask(mask=('[#3e ]', ), )
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
' -------------------------------  Assembly -------------------------------------' 
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
' -------------------------------  Step -------------------------------------' 
# *****************************************************************************    
    
def Step():
    """ABAQUS_function: Step:
           Return: Stress, Strain, Elastic energy """
   
           
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
            'S', 'E', 'U', 'IVOL', 'ELSE'))

    
# =============================================================================
# F-Functions    
# =============================================================================
    
def F1_Geometry(c,t,h):
    '''F1_Geometry: Designs model based on 4 independents geometrical parameters; c,t,h
        c: length of the internal hexagonal.
        t: the thickness of the honeycomb.
        tg: tg = t it's just to avoid errors in a function of Fill_hex_external.
        h: the hight of the honeycomb.
        X, Y, Z: Cube size.
        pi: points of internal hexagonal. (i \in {1,...,6})
        qi: points of external hexagonal. (i \in {1,...,4})
    '''
    #Calculation other Parameters required-------------------------------------------------------------
    (X,Z,Y) = f_cth_to_XZY(c,t,h)
    (p1, p2, p3, p4, p5, p6) = hexagonale_points(c,t,X,Z)               # To cut the internal hexagonal
    (q1, q2, q3, q4)= hexagonale_points_External(c,t,X,Z)               # To cut the external hexagonal 
    #Model----------------------------------------------------------------------------------------------
    Cube(X,Y,Z)                                                         # Make a cube with size X, Y, Z
    Cut_hex(h, Y, Z, p1, p2, p3, p4, p5, p6)                           # Cut the internal Hexagonal 
    Cut_hex_external(h, Y, Z, q1, q2, q3, q4)                          # Cut the external Hexagonal           
    Fill_hex(h, Y, Z, t, p1, p2, p3, p4, p5, p6)                       # Fill the internal Hexagonal with air!
    Fill_hex_external(h, Y, Z, t, q1, q2, q3, q4)                     # Fill the external Hexagonal with air!

def F2_Materials(E_Solid, v_Solid, name_Solid, E_Air, v_Air, name_Air):
    '''F2_Materials: Considers elastic modules for solid and air (with their names)
        E_Solid, E_Air: Elastic modules of solid and void parts.
        v_Solid, v_Solid: poisson ratio of solid and void parts
        name_Solid, name_Air: names of materials 
    '''    
    Material_Section_Solid(E_Solid, v_Solid, name_Solid, name_Solid)   # Solid Part
    Material_Section_Void(E_Air, v_Air, name_Air, name_Air)            # Void Part
    Material_color()   

def F3_Assembly(X,Y,Z):
    'F3_Assembly: Assemble and move the coordinate to the center of model'
    Assembly()                 # Assembly the model
    Move_Ref_to(0.0, 0.0, -Z)  # Gives the origin of coordinates (i.e. 0,0,0) at the center of model
# =============================================================================
# =============================================================================
# =============================================================================
# # # MESH 
# =============================================================================
# =============================================================================
# =============================================================================
def edge_seed_fb(edge_fb, dF):

    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['Part-1'].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#14018030 #a0 #303 ]', ), )
    a.seedEdgeBySize(edges=pickedEdges, size=edge_fb, deviationFactor=dF, 
        constraint=FINER)


def edge_seed_rl(edge_rl, dF):

    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['Part-1'].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#1402805 #f000000a ]', ), )
    a.seedEdgeBySize(edges=pickedEdges, size=edge_rl, deviationFactor=dF, 
        constraint=FINER)


def edge_seed_tb(edge_tb, dF):

    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['Part-1'].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#0:2 #821040 ]', ), )
    a.seedEdgeBySize(edges=pickedEdges, size=edge_tb, deviationFactor=dF,
                     constraint=FINER)


def partition_cell():
    'Defines partition for havin better mesh'

   
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#1 ]', ), )
    v11 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v11[15], point2=v11[31], 
        point3=v11[30], cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#2 ]', ), )
    v1 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v1[0], point2=v1[41], point3=v1[7], 
        cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#4 ]', ), )
    v11 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v11[17], point2=v11[3], point3=v11[0], 
        cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#2 ]', ), )
    v1 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v1[25], point2=v1[17], point3=v1[26], 
        cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#8 ]', ), )
    v11 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v11[5], point2=v11[26], 
        point3=v11[17], cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#10 ]', ), )
    v1 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v1[0], point2=v1[13], point3=v1[4], 
        cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#2 ]', ), )
    v11 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v11[23], point2=v11[13], 
        point3=v11[22], cells=pickedCells)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedCells = c1.getSequenceFromMask(mask=('[#20 ]', ), )
    v1 = a.instances['Part-1'].vertices
    a.PartitionCellByPlaneThreePoints(point1=v1[30], point2=v1[28], point3=v1[31], 
        cells=pickedCells)



def Mesh_HEX_DOMINATED_Q_Advance_front(element_size,min_size_fac, dF):

    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    pickedRegions = c1.getSequenceFromMask(mask=('[#7fff ]', ), )
    a.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED, 
        technique=SWEEP, algorithm=ADVANCING_FRONT)
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Part-1'].cells
    cells1 = c1.getSequenceFromMask(mask=('[#7fff ]', ), )
    pickedRegions =(cells1, )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))
    
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.seedPartInstance(regions=partInstances, size=element_size, deviationFactor= dF, 
        minSizeFactor=min_size_fac)
    
    a = mdb.models['Model-1'].rootAssembly
    partInstances =(a.instances['Part-1'], )
    a.generateMesh(regions=partInstances)








