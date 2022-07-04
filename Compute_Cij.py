# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:03:30 2019

@author: momoe

This code computes the effective properties of an RVE using numerical homogenization
"""

import numpy as np
import time

# *****************************************************************************
' -------------------------------  Imports ------------------------------------' 
# ***************************************************************************** 
' Abaqus modules '
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


# *****************************************************************************
' -------------------------------  Sets -------------------------------------' 
# *****************************************************************************  
# =============================================================================
# Walls
# =============================================================================
def Functions_Set_Face_AC(Face, Name_Face):
    """
    Assigns a face.
    (To reduce the lines and volume of each set function)
    
    Input:
       Face: set of nodes (list) 
       Name_Face: name of that set of nodes (string)
    """   
         
    modelName1='Model-1'
    mdb.models[modelName1].rootAssembly.Set(faces=Face, name=Name_Face)
    mdb.models[modelName1].rootAssembly.Set(name=Name_Face,
              nodes=mdb.models[modelName1].rootAssembly.sets[Name_Face].nodes)

    

def Functions_Set_walls(X,Y,Z):
     """
     Assigns 6 sets of nodes for 6 surfaces of the RVE, namely as:
         Top/Bottom, Front/Back, Right/Left
         
    Input:
        X, Y, Z: coordinate of the one corner of the cube (positive scalar)
        
                                    --------
                                   /       /|
                                  /       / |
                                 ---------  |
                                 |       |  |
                                 |       | / 
                                 ---------/
        
    """     
         
     modelName1='Model-1'
     instanceName1='Part-1'      
    #The values should be round-off
     X=round(X,6); Y=round(Y,6); Z=round(Z,6);
      
     FRONT=[]; BACK=[];
     LEFT=[]; RIGHT=[];
     TOP=[]; BOTTOM=[];
     ALL_SURFACE_NODES=[]
    

     for i in mdb.models[modelName1].rootAssembly.instances[instanceName1].faces:
         a=i.pointOn[0]    
    
         if a[2]== Z:
            FRONT=FRONT+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]
         elif a[2]== -Z:
             BACK=BACK+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]
         elif a[0]== X:
             RIGHT=RIGHT+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]
         elif a[0]== -X:
        	 LEFT=LEFT+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]
         elif a[1]== -Y:
        	 BOTTOM=BOTTOM+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]        
         elif a[1]== Y:
        	 TOP=TOP+[mdb.models[modelName1].rootAssembly.instances[instanceName1].faces.findAt(((a[0],a[1],a[2]),))]

    #Assign set
     Functions_Set_Face_AC(FRONT, 'FRONT')
     Functions_Set_Face_AC(BACK, 'BACK')
     Functions_Set_Face_AC(LEFT, 'LEFT')
     Functions_Set_Face_AC(RIGHT, 'RIGHT')
     Functions_Set_Face_AC(TOP, 'TOP')
     Functions_Set_Face_AC(BOTTOM, 'BOTTOM')



# =============================================================================
# Vertices
# =============================================================================
def Functions_Set_vertices(X,Y,Z, tole_Ver = 0.001, 
                           ac_x = 4, ac_y = 3, ac_z = 4):
    """
    Assigns a set for each vertices of the cube.  
    
    Input:
        X, Y, Z: coordinate of the one corner of the cube (positive scalar)
        tole_Ver: tolerance to find the vertices. 
        ac_x, ac_y, ac_z: numerical precision to tound the coordinate values
    
    Note:
        For complex geometry (eg.: low relative density in cellular materials), 
        probably the input values tole_Ver, ac_x, ac_y and ac_z should be changed to avoid
        the error. Otherwise, this fucntion cannot find the vertices and once
        the next function wants to put a kinematic constraint (PBC), it will
        causes an error. 
        
    """
    
    modelName1='Model-1'
    instanceName1='Part-1'
    all_vertices=[]
    X= round(X, ac_x); 
    Y= round(Y, ac_y); 
    Z= round(Z, ac_z);
    
    def Set_point(node_X,node_Y,node_Z, X0, Y0, Z0, Name_vertice):
        modelName1='Model-1'
        instanceName1='Part-1'
        point_info=[]
        if abs(node_X-X0)<tole_Ver and abs(node_Y-Y0)<tole_Ver and abs(node_Z-Z0)<tole_Ver:
            
            point_info=point_info+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
            mdb.models[modelName1].rootAssembly.Set(name= Name_vertice, nodes=
                      mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(point_info[0][3]-1):(point_info[0][3])])
 
            
    for i in mdb.models[modelName1].rootAssembly.instances[instanceName1].nodes:
        node_X=round(i.coordinates[0],4); 
        node_Y=round(i.coordinates[1],3);
        node_Z=round(i.coordinates[2],4);
        
        
        try:
            Set_point(node_X,node_Y,node_Z, X, -Y, Z, 'V_FrRiBo')
            Set_point(node_X,node_Y,node_Z, -X, -Y, Z, 'V_FrLeBo')
            Set_point(node_X,node_Y,node_Z, X, +Y, Z, 'V_FrRiTo')
            Set_point(node_X,node_Y,node_Z, -X, +Y, Z, 'V_FrLeTo')
                
            Set_point(node_X,node_Y,node_Z, -X, -Y, -Z, 'V_BaLeBo')
            Set_point(node_X,node_Y,node_Z, X, -Y, -Z, 'V_BaRiBo')
            Set_point(node_X,node_Y,node_Z, -X, +Y, -Z, 'V_BaLeTo')
            Set_point(node_X,node_Y,node_Z, X, +Y, -Z, 'V_BaRiTo')
            
        except:
            print("The Functions_Set_vertices could not assign the vertices.")
            raise ValueError
  

        
# =============================================================================
# Edges
# =============================================================================
def Functions_Set_Edges(face_1, face_2, edge_name, X,Y,Z,
                        ac_x = 4, ac_y = 3, ac_z = 4): 
    """
    Assigns a set for each eades of the cube.  
    
    Input:
        X, Y, Z: coordinate of the one corner of the cube (positive scalar)
        edge_name: name of the edge (string)
        face_1, face_2: mutual faces
        tole_Ver: tolerance to find the vertices. 
        ac_x, ac_y, ac_z: numerical precision to tound the coordinate values
    
    Note:
        This function is working with the node's label. 
        For complex geometry (low relative density in my case), yet, it causes 
        an error. The input values tole_Ver, ac_x, ac_y and ac_z should be 
        changed to avoid the error. Otherwise, this fucntion cannot find the
        vertices and once the next function wants to put a kinematic
        constraint (PBC), it will cause an error. 
        
    """
    modelName1='Model-1'
    instanceName1='Part-1'
    X= round(X, ac_x); Y= round(Y, ac_y); Z= round(Z, ac_z); 
    c=[]
    for i in mdb.models[modelName1].rootAssembly.sets[face_1].nodes:
        node_X=abs(round(i.coordinates[0], ac_x));
        node_Y=abs(round(i.coordinates[1], ac_y)); 
        node_Z=abs(round(i.coordinates[2], ac_z));
        # abs() to check all vertices, round to control precision
        for j in mdb.models[modelName1].rootAssembly.sets[face_2].nodes:
            if i.label == j.label and [node_X, node_Y, node_Z] != [X,Y,Z]:
                c=c+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
                                     
    c.sort()
    rep=1
    for i in c:
        mdb.models[modelName1].rootAssembly.Set(name= edge_name+ '-' +str(rep), nodes=
                 mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
        rep=rep+1
    
    number_of_edge_nodes=rep-1
    'You would need number number_of_edge_nodes for constraint edges loop' 
    
    return number_of_edge_nodes, c    

# =============================================================================
# dummy1,2,3 
# =============================================================================
def Set_dummy_points(X,Y,Z, cx,cy,cz, Name):
    """
    Assigns a three dummy (reference) points to apply the displacement. 
    
    Input:
        X, Y, Z: coordinate of the one corner of the cube (positive scalar)
        cx,cy,cz: coefficients to provide the dummy point coordinate:
            (x, y, z)_dummy = (X*cx,  Y*cy, Z*cz)
    
    Note:
        Input the (cx, cy, cz) > 1.0 then the dummy points will be beside the 
        surfaces. 
    
    """
    mdb.models['Model-1'].rootAssembly.AttachmentPoints(name= Name+'1', 
              points=((X*cx, 0.0, 0.0), ), setName= Name+'1')

    mdb.models['Model-1'].rootAssembly.AttachmentPoints(name= Name+'2', 
              points=((0.0, Y*cy, 0.0), ), setName= Name+'2')

    mdb.models['Model-1'].rootAssembly.AttachmentPoints(name= Name+'3', 
              points=((0.0, 0.0, Z*cz), ), setName= Name+'3')
    
    if cx < 1.0 or cy < 1.0 or cz < 1.0:
        print('*** Note: It is better to have (cx, cy, cz) > 1.0 ***')
  
# =============================================================================
# Edge
# =============================================================================
def face_without_edges(Face_name, all_edges):
    """
    Return bx which includes nodes of the face without edges
    
    Input:
       Face_name: name of the set of nodes over the walls (Front, Top, ...) 
           
       all_edges: a list includes the nodes over the all edges of the cube.
      
    Note:
        This function has not been used. 
    """ 
    modelName1='Model-1'
    b=[]
    for i in mdb.models[modelName1].rootAssembly.sets[Face_name].nodes:
        b=b+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b.sort()

    bx=[]
    for i in b:
        if i not in all_edges:
            bx.append(i)                                   
    bx.sort()
    
    return bx

# =============================================================================
# closest point
# =============================================================================
def Select_closet_node_to_d0(d0x,d0y,d0z,set_name, tole_centerP):
    """
    Assigns the closest point to a given specific node d0(d0x, d0y, d0z)
    with telorance telo.
    
    d0x,d0y,d0z: coordinate of the target node
    set_name: name of the set.
    tole_centerP: tolerance 
    
    """
    
    modelName1='Model-1' 
    instanceName1='Part-1'
    center_point = []

    for i in mdb.models[modelName1].rootAssembly.instances[instanceName1].nodes:
        d=(i.coordinates[0]-d0x)**2 + (i.coordinates[1]-d0y)**2 + (i.coordinates[2]-d0z)**2
        if d < tole_centerP:
            center_point=[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
            
        elif tole_centerP <d < 4.*tole_centerP:
            center_point_2=[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    
    if len(center_point) == 0:
        print('Warning: The tolerance tole_centerP was too smal. It is increased to 4.*tole_centerP = ', 4.*tole_centerP)
        for i in center_point_2:
                mdb.models[modelName1].rootAssembly.Set(name=set_name, nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                
    else:
        for i in center_point:
                mdb.models[modelName1].rootAssembly.Set(name=set_name, nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
    
    
            
        

# *****************************************************************************
' ----------------------  Constraints functions -------------------------------' 
# *****************************************************************************

# =============================================================================
# AXIAL
# =============================================================================
def Constraint_RIGHT_LEFT_Axial(element_size, telo, X,Y,Z):
    '''
    This function defines two lists including the coordinates and label of the
    nodes over the LEFT & RIGHT faces. Then, put constraint between two
    corresponding nodes with one condition;
        if the distance of two nodes is less than the "telo*element_size"
    
    In this way, however, it is possible that for some nodes the fucntion can't
    find the corresponding nodes or, even worse, find more than one node for another
    nodes. This casues the "missing degree of freedom" error and the analysis 
    will be aborted. In generall three cases can be happened:
    
    "OVER-CONSTRAINED": When one node has constraints with two (or more) nodes.
    "Perfect-CONSTRAINED": When one node has constraints with one (and only one) other node.
    "UNDER-CONSTRAINED": When one node has no constraints.
    
    Input:
        X, Y, Z: coordinate of the one corner of the RVE (positive scalar)
        element_size: element size
        telo: tolerance to find the correspondig node. 
        
    Parameters:
        previous_node: I defined this parameters to avoid repetetive constraint
        (over_constrained)
    
    Note and tips:
        Once you have some very distorted element compare to other element, 
        it wiuld be hard to guess the 'telo'. It is hard becasue if you put it
        high, it will be OVER-CONSTRAINED (find more than one) and if you put
        it low it would be 'UNDER-CONSTRAINED' (it can't find any node). In my
        RVE, it happens when the relative density is too low (or c/t is too high).
        
        To avoid over constrain I wrote 'i != previous_node' in the for loop. 
        
    
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_LEFT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['LEFT'].nodes:
        if round(Z,4)!=abs(round(i.coordinates[2],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):
            b_LEFT=b_LEFT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_LEFT.sort()
    
    b_RIGHT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['RIGHT'].nodes:
        if round(Z,4)!=abs(round(i.coordinates[2],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):

            b_RIGHT=b_RIGHT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_RIGHT.sort()
    
    nodes_number_right = len(b_RIGHT)
    nodes_number_left = len(b_LEFT)
    active_nodes_RL = 0
        
    rep=1
    x=0
    previous_node='not the first one'
    previous_node_left='not the first one'
    selected_nodes_R = []
    selected_nodes_L = []
    # i = [X, Y, Z, label]
    # We don't care about the first point. it will start from second one
    for i in b_RIGHT:
        for j in b_LEFT:
#            if abs(i[1]-j[1])<telo*element_size and abs(i[2]-j[2])<telo*element_size and i != previous_node and j != previous_node_left:
            if abs(i[1]-j[1])<telo*element_size and abs(i[2]-j[2])<telo*element_size and i[3] not in selected_nodes_R and j[3] not in selected_nodes_L:
                # Part_I and Part_II: Geometrical conditions
                # Part_III: To avoid overconstrained 
                active_nodes_RL += 1
                
                mdb.models[modelName1].rootAssembly.Set(name='CANode-RIGHT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-LEFT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-R&L-x-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 1),(1.0,'CANode-LEFT-'+str(rep), 1), 
                                 (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-R&L-y-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 2),(1.0,'CANode-LEFT-'+str(rep), 2), 
                                 (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-R&L-z-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 3),(1.0,'CANode-LEFT-'+str(rep), 3), 
                                 (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint
                previous_node_left = j #We don't want to put repetitive constraint
                selected_nodes_R.append(i[3])
                selected_nodes_L.append(j[3])
                
            
                

    print('nodes_number_right=', nodes_number_right)
    print('nodes_number_left=', nodes_number_left)
    print('active_nodes_RL=', active_nodes_RL)
    
    # uncomment the following line for complex RVE
#    if nodes_number_right < active_nodes_RL or nodes_number_left < active_nodes_RL:
#        print('****** WARNING: OVER-CONSTRAINED ---> Left and Right******')
#
#    if nodes_number_right > active_nodes_RL or nodes_number_left > active_nodes_RL:
#        print('****** WARNING: UNDER-CONSTRAINED ---> Left and Right ******')
#        print('There are: ', nodes_number_right - active_nodes_RL , 'nodes WITHOUT constraints')
#        raise ValueError
        
    if nodes_number_right == active_nodes_RL or nodes_number_left == active_nodes_RL:
        print('____ Perfect CONSTRAINED ---> Left and Right _____')
        



def Constraint_FRONT_BACK_Axial(element_size, telo, X,Y,Z):
    '''
    This function firs defins two sets to get coordinate and label of each node of FRONT & BACK faces. 
    Then, put constraint between two corresponding nodes with one condition; if the distance of two nodes is 
    less than telo*element_size
    
    OVER-CONSTRAINED: When one node has constraints with two (or more) nodes.
    Perfect-CONSTRAINED: When one node has constraints with one (and only one) other node.
    UNDER-CONSTRAINED: When one node has no constraints.
    
    previous_node: I defined this parameters to avoid repetetive constraint (over_constrained)
    
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_FRONT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['FRONT'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):
            b_FRONT=b_FRONT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_FRONT.sort()
    
    b_BACK=[]
    for i in mdb.models[modelName1].rootAssembly.sets['BACK'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):

            b_BACK=b_BACK+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_BACK.sort()
    
    nodes_number_front = len(b_FRONT)
    nodes_number_back = len(b_BACK)
    active_nodes_FB = 0
        
    rep=1
    x=0
    previous_node='not the first one' # We don't care about the first point. it will start from second one
    selected_nodes_F = []
    selected_nodes_B = []
    for i in b_FRONT:
        for j in b_BACK:
#            if abs(i[0]-j[0])<telo*element_size and abs(i[1]-j[1])<telo*element_size  and i != previous_node:
            if abs(i[0]-j[0])<telo*element_size and abs(i[1]-j[1])<telo*element_size  and i[3] not in selected_nodes_F and j[3] not in selected_nodes_B:
                # Part_I and Part_II: Geometrical condition ----  Part_III: To avoid overconstrained
                active_nodes_FB +=  1
                mdb.models[modelName1].rootAssembly.Set(name='CANode-FRONT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-BACK-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-F&B-x-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 1),(1.0,'CANode-BACK-'+str(rep), 1), 
                                 (0.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-F&B-y-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 2),(1.0,'CANode-BACK-'+str(rep), 2), 
                                 (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-F&B-z-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 3),(1.0,'CANode-BACK-'+str(rep), 3), 
                                 (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint
                selected_nodes_F.append(i[3])
                selected_nodes_B.append(j[3])
                
                '''Note: It is true that we do nor have periodic boundary condition in Y-direction,
                But If we delet the Y-constraint we woud have a distorsion in sides (specially in front and back).
                Try it if you wanna see!
                '''

        
    print('nodes_number_front=', nodes_number_front)
    print('nodes_number_back=', nodes_number_back)
    print('active_nodes_FB=', active_nodes_FB)
    
    if nodes_number_front < active_nodes_FB or nodes_number_back < active_nodes_FB:
        print('****** WARNING: OVER-CONSTRAINED ---> Front and Back******')

    if nodes_number_front > active_nodes_FB or nodes_number_back > active_nodes_FB:

        print('****** WARNING: UNDER-CONSTRAINED ---> Front and Back ******')
        print('There are: ', nodes_number_front - active_nodes_FB , 'nodes WITHOUT constraints')
        
    if nodes_number_front == active_nodes_FB or nodes_number_back == active_nodes_FB:
        print('____ Perfect CONSTRAINED ---> Front and Back _____')




def Constraint_TOP_BOTTOM_Axial(element_size, telo, X,Y,Z):
    '''
    This function firs defins two sets to get coordinate and label of each node of FRONT & BACK faces. 
    Then, put constraint between two corresponding nodes with one condition; if the distance of two nodes is 
    less than telo*element_size
    
    OVER-CONSTRAINED: When one node has constraints with two (or more) nodes.
    Perfect-CONSTRAINED: When one node has constraints with one (and only one) other node.
    UNDER-CONSTRAINED: When one node has no constraints.
    
    previous_node: I defined this parameters to avoid repetetive constraint (over_constrained)
    
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_TOP=[]
    for i in mdb.models[modelName1].rootAssembly.sets['TOP'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Z,4)!=abs(round(i.coordinates[2],4)):
            b_TOP=b_TOP+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_TOP.sort()
    
    b_BOTTOM=[]
    for i in mdb.models[modelName1].rootAssembly.sets['BOTTOM'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Z,4)!=abs(round(i.coordinates[2],4)):

            b_BOTTOM=b_BOTTOM+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_BOTTOM.sort()
    
    nodes_number_top = len(b_TOP)
    nodes_number_bottom = len(b_BOTTOM)
    active_nodes_TB = 0
        
    rep=1
    x=0
    previous_node='not the first one' # We don't care about the first point. it will start from second one
    selected_nodes_T = []
    selected_nodes_BO = []
    for i in b_TOP:
        for j in b_BOTTOM:
            if abs(i[0]-j[0])<telo*element_size and abs(i[2]-j[2])<telo*element_size  and i[3] not in selected_nodes_T and j[3] not in selected_nodes_BO:
                # Part_I and Part_II: Geometrical condition ----  Part_III: To avoid overconstrained
                active_nodes_TB +=  1
                mdb.models[modelName1].rootAssembly.Set(name='CANode-TOP-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-BOTTOM-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-T&B-x-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 1),(1.0,'CANode-BOTTOM-'+str(rep), 1), 
                                 (0.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-T&B-y-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 2),(1.0,'CANode-BOTTOM-'+str(rep), 2), 
                                 (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-T&B-z-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 3),(1.0,'CANode-BOTTOM-'+str(rep), 3), 
                                 (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint
                selected_nodes_T.append(i[3])
                selected_nodes_BO.append(j[3])
                
                '''Note: It is true that we do nor have periodic boundary condition in Y-direction,
                But If we delet the Y-constraint we woud have a distorsion in sides (specially in front and back).
                Try it if you wanna see!
                '''

        
    print('nodes_number_top=', nodes_number_top)
    print('nodes_number_bottom=', nodes_number_bottom)
    print('active_nodes_TB=', active_nodes_TB)
    
    if nodes_number_top < active_nodes_TB or nodes_number_bottom < active_nodes_TB:
        print('****** WARNING: OVER-CONSTRAINED ---> Top and Bottom******')

    if nodes_number_top > active_nodes_TB or nodes_number_bottom > active_nodes_TB:

        print('****** WARNING: UNDER-CONSTRAINED ---> Top and Bottom ******')
        print('There are: ', nodes_number_top - active_nodes_TB , 'nodes WITHOUT constraints')
        
    if nodes_number_top == active_nodes_TB or nodes_number_bottom == active_nodes_TB:
        print('____ Perfect CONSTRAINED ---> Top and Bottom _____')



# =============================================================================
# # Egdes - Constraints AXIAL: E11, E22, E33 (to avoid overlaping)
# =============================================================================
def ConstraintsEdgesAxial_FrRi_BaLe(last_node_FRONT_LEFT):
    modelName1='Model-1'
    "Apply MPC over edges: E_FrRi=E_BaLe"  
    N=last_node_FRONT_LEFT        
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_FrRi=E_BaLe}-'+str(i+1), 
                 terms=((-1.0,'E_FrRi-'+str(rep), 1),(1.0,'E_BaLe-'+str(rep), 1),
                        (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_FrRi=E_BaLe}'+str(i+1), 
    	    terms=((-1.0,'E_FrRi-'+str(rep), 2),(+1.0,'E_BaLe-'+str(rep), 2), 
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1   
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_FrRi=E_BaLe}'+str(i+1), 
                  terms=((-1.0,'E_FrRi-'+str(rep), 3),(1.0,'E_BaLe-'+str(rep), 3),
                    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
        rep=rep+1

#------------------------------------------------------------------------------
def ConstraintsEdgesAxial_BaRi_FrLe(last_node_BACK_RIGHT):       
    "Apply MPC over edges: E_BaRi=E_FrLe"
    modelName1='Model-1'
    N=last_node_BACK_RIGHT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaRi=E_FrLe}-'+str(i+1), 
                  terms=((-1.0,'E_BaRi-'+str(rep), 1),(1.0,'E_FrLe-'+str(rep), 1),
                    (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
                    
        rep=rep+1
    #-------    
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaRi=E_FrLe}'+str(i+1), 
    	    terms=((-1.0,'E_BaRi-'+str(rep), 2),(1.0,'E_FrLe-'+str(rep), 2), 
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1 
    #-------             
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaRi=E_FrLe}'+str(i+1), 
            terms=((1.0,'E_BaRi-'+str(rep), 3),(-1.0,'E_FrLe-'+str(rep), 3),
                   (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
        rep=rep+1 

#------------------------------------------------------------------------------
def ConstraintsEdgesAxial_ToLe_BoRi(last_node_TOP_LEFT):
    "Apply MPC over edges: E_ToLe=E_BoRi"
    modelName1='Model-1'
    N=last_node_TOP_LEFT    
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((1.0,'E_ToLe-'+str(rep), 1),(-1.0,'E_BoRi-'+str(rep), 1),
                (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------  
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
    	    terms=((-1.0,'E_ToLe-'+str(rep), 2), (1.0,'E_BoRi-'+str(rep), 2),
	       (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))
        rep=rep+1
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((-1.0,'E_ToLe-'+str(rep), 3),(1.0,'E_BoRi-'+str(rep), 3),
	        (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1

#------------------------------------------------------------------------------
def ConstraintsEdgesAxial_BoLe_ToRi(last_node_TOP_LEFT):          
    "Apply MPC over edges: E_BoLe=E_ToRi"
    modelName1='Model-1'
    N=last_node_TOP_LEFT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
		    terms=((1.0,'E_BoLe-'+str(rep), 1),(-1.0,'E_ToRi-'+str(rep), 1),
     		(1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------          
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
    	    terms=((1.0,'E_BoLe-'+str(rep), 2),(-1.0,'E_ToRi-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1    
    #-------   
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
     	    terms=((-1.0,'E_BoLe-'+str(rep), 3),(1.0,'E_ToRi-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
    	rep=rep+1    

#------------------------------------------------------------------------------
def ConstraintsEdgesAxial_ToFr_BaBo(last_node_TOP_FRONT):
    "Apply MPC over edges:  E_ToFr=E_BaBo"
    modelName1='Model-1'       
    N=last_node_TOP_FRONT 
    rep=1
    for i in range(0,N):	
    	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	    	terms=((-1.0,'E_ToFr-'+str(rep), 1),(1.0,'E_BaBo-'+str(rep), 1),
		    (0.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------   
    rep=1
    for i in range(0,N):
	    mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
    	    terms=((-1.0,'E_ToFr-'+str(rep), 2), (1.0,'E_BaBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2))); rep=rep+1
        
    #-------
    rep=1   
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	       terms=((-1.0,'E_ToFr-'+str(rep), 3),(1.0,'E_BaBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
        rep=rep+1

#------------------------------------------------------------------------------
def ConstraintsEdgesAxial_BaTo_FrBo(last_node_TOP_FRONT):
    "Apply MPC over edges: E_BaTo=E_FrBo"
    modelName1='Model-1'
    N=last_node_TOP_FRONT
    rep=1
    for i in range(0,N):	
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    		terms=((1.0,'E_BaTo-'+str(rep), 1),(-1.0,'E_FrBo-'+str(rep), 1),
    		(0.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 2),(1.0,'E_FrBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(1.0,'dummyV2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1
    rep=1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((1.0,'E_BaTo-'+str(rep), 3),(-1.0,'E_FrBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
    	rep=rep+1

# =============================================================================
# # VERTICES - Constraints AXIAL (to avoid overlaping)
# =============================================================================

def Constraints_Vertices_Axial():
    modelName1='Model-1' 
    #X-direction - DOF=1
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((1.0, 'V_BaLeBo', 1),(-1.0, 'V_FrRiTo', 1), (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((1.0, 'V_BaLeTo', 1),(-1.0, 'V_FrRiBo', 1), (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((1.0, 'V_FrLeTo', 1),(-1.0, 'V_BaRiBo', 1), (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 1),(-1.0, 'V_BaRiTo', 1), (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))

    #Y-direction - DOF=2-----------------------------------------------------------
    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((1.0, 'V_BaLeBo', 2),(-1.0, 'V_FrRiTo', 2), (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 2),(1.0, 'V_FrRiBo', 2), (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 2),(1.0, 'V_BaRiBo', 2), (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 2),(-1.0, 'V_BaRiTo' , 2), (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))

    #Z-direction - DOF=3-----------------------------------------------------------
    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((1.0, 'V_BaLeBo', 3),(-1.0, 'V_FrRiTo', 3), (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((1.0, 'V_BaLeTo', 3),(-1.0, 'V_FrRiBo', 3), (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-3}', 
     	terms=((1.0, 'V_BaRiBo', 3),(-1.0, 'V_FrLeTo' , 3), (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_BaRiTo' , 3),(-1.0, 'V_FrLeBo', 3), (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))



# =============================================================================
# =============================================================================
# # Shear
# =============================================================================
# =============================================================================
def Constraint_RIGHT_LEFT_Shear(element_size, telo, X,Y,Z):
    '''
    This function first defins two sets to get coordinate and label of each node of LEFT & RIGHT faces. 
    Then, put constraint between two corresponding nodes with one condition; if the distance of two nodes is 
    less than the "telo*element_size"
    
    "OVER-CONSTRAINED": When one node has constraints with two (or more) nodes.
    "Perfect-CONSTRAINED": When one node has constraints with one (and only one) other node.
    "UNDER-CONSTRAINED": When one node has no constraints.
    
    previous_node: I defined this parameters to avoid repetetive constraint (over_constrained)
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_LEFT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['LEFT'].nodes:
        if round(Z,4)!=abs(round(i.coordinates[2],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):
            b_LEFT=b_LEFT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_LEFT.sort()
    
    b_RIGHT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['RIGHT'].nodes:
        if round(Z,4)!=abs(round(i.coordinates[2],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):

            b_RIGHT=b_RIGHT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_RIGHT.sort()
    
    nodes_number_right = len(b_RIGHT)
    nodes_number_left = len(b_LEFT)
    active_nodes_RL = 0
        
    rep=1
    x=0
    previous_node='not the first one' # We don't care about the first point. it will start from second one
    for i in b_RIGHT:
        for j in b_LEFT:
            if abs(i[1]-j[1])<telo*element_size and abs(i[2]-j[2])<telo*element_size and i != previous_node:
                # Part_I and Part_II: Geometrical condition ----  Part_III: To avoid overconstrained 
                active_nodes_RL += 1
                
                mdb.models[modelName1].rootAssembly.Set(name='CANode-RIGHT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-LEFT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-R&L-x-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 1),(1.0,'CANode-LEFT-'+str(rep), 1), 
                                 (1.0,'dummy1', 1),(0.0,'dummy2', 1),(0.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-R&L-y-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 2),(1.0,'CANode-LEFT-'+str(rep), 2), 
                                 (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-R&L-z-'+str(x+1), 
                          terms=((-1.0,'CANode-RIGHT-'+str(rep), 3),(1.0,'CANode-LEFT-'+str(rep), 3), 
                                 (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint 
                
                '''Note: If we delete the Y-constraint we woud have a distorsion in sides (specially in front and back).
                Try it if you wanna see!
                '''


    print('nodes_number_right=', nodes_number_right)
    print('nodes_number_left=', nodes_number_left)
    print('active_nodes_RL=', active_nodes_RL)
    
    if nodes_number_right < active_nodes_RL or nodes_number_left < active_nodes_RL:
        print('****** WARNING: OVER-CONSTRAINED ---> Left and Right******')

    if nodes_number_right > active_nodes_RL or nodes_number_left > active_nodes_RL:
        print('****** WARNING: UNDER-CONSTRAINED ---> Left and Right ******')
        print('There are: ', nodes_number_right - active_nodes_RL , 'nodes WITHOUT constraints')
        
    if nodes_number_right == active_nodes_RL or nodes_number_left == active_nodes_RL:
        print('____ Perfect CONSTRAINED ---> Left and Right _____')
        


def Constraint_FRONT_BACK_Shear(element_size, telo, X,Y,Z):
    '''
    This function firs defins two sets to get coordinate and label of each node of FRONT & BACK faces. 
    Then, put constraint between two corresponding nodes with one condition; if the distance of two nodes is 
    less than telo*element_size
    
    OVER-CONSTRAINED: When one node has constraints with two (or more) nodes.
    Perfect-CONSTRAINED: When one node has constraints with one (and only one) other node.
    UNDER-CONSTRAINED: When one node has no constraints.
    
    previous_node: I defined this parameters to avoid repetetive constraint (over_constrained)
    
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_FRONT=[]
    for i in mdb.models[modelName1].rootAssembly.sets['FRONT'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):
            b_FRONT=b_FRONT+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_FRONT.sort()
    
    b_BACK=[]
    for i in mdb.models[modelName1].rootAssembly.sets['BACK'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Y,4)!=abs(round(i.coordinates[1],4)):

            b_BACK=b_BACK+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_BACK.sort()
    
    nodes_number_front = len(b_FRONT)
    nodes_number_back = len(b_BACK)
    active_nodes_FB = 0
        
    rep=1
    x=0
    previous_node='not the first one' # We don't care about the first point. it will start from second one
    for i in b_FRONT:
        for j in b_BACK:
            if abs(i[0]-j[0])<telo*element_size and abs(i[1]-j[1])<telo*element_size  and i != previous_node:
                # Part_I and Part_II: Geometrical condition ----  Part_III: To avoid overconstrained
                active_nodes_FB +=  1
                mdb.models[modelName1].rootAssembly.Set(name='CANode-FRONT-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-BACK-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-F&B-x-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 1),(1.0,'CANode-BACK-'+str(rep), 1), 
                                 (0.0,'dummy1', 1),(0.0,'dummy2', 1),(1.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-F&B-y-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 2),(1.0,'CANode-BACK-'+str(rep), 2), 
                                 (0.0,'dummy1', 2),(0.0,'dummy2', 2),(1.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-F&B-z-'+str(x+1), 
                          terms=((-1.0,'CANode-FRONT-'+str(rep), 3),(1.0,'CANode-BACK-'+str(rep), 3), 
                                 (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint
                
                '''Note: If we delet the Y-constraint we woud have a distorsion in sides (specially in front and back).
                Try it if you wanna see!
                '''

        
    print('nodes_number_front=', nodes_number_front)
    print('nodes_number_back=', nodes_number_back)
    print('active_nodes_FB=', active_nodes_FB)
    
    if nodes_number_front < active_nodes_FB or nodes_number_back < active_nodes_FB:
        print('****** WARNING: OVER-CONSTRAINED ---> Front and Back******')

    if nodes_number_front > active_nodes_FB or nodes_number_back > active_nodes_FB:

        print('****** WARNING: UNDER-CONSTRAINED ---> Front and Back ******')
        print('There are: ', nodes_number_front - active_nodes_FB , 'nodes WITHOUT constraints')
        
    if nodes_number_front == active_nodes_FB or nodes_number_back == active_nodes_FB:
        print('____ Perfect CONSTRAINED ---> Front and Back _____')




def Constraint_TOP_BOTTOM_Shear(element_size, telo, X,Y,Z):
    '''
    This function firs defins two sets to get coordinate and label of each node of FRONT & BACK faces. 
    Then, put constraint between two corresponding nodes with one condition; if the distance of two nodes is 
    less than telo*element_size
    
    OVER-CONSTRAINED: When one node has constraints with two (or more) nodes.
    Perfect-CONSTRAINED: When one node has constraints with one (and only one) other node.
    UNDER-CONSTRAINED: When one node has no constraints.
    
    previous_node: I defined this parameters to avoid repetetive constraint (over_constrained)
    
    '''
    modelName1='Model-1'
    instanceName1='Part-1'
    
    b_TOP=[]
    for i in mdb.models[modelName1].rootAssembly.sets['TOP'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Z,4)!=abs(round(i.coordinates[2],4)):
            b_TOP=b_TOP+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_TOP.sort()
    
    b_BOTTOM=[]
    for i in mdb.models[modelName1].rootAssembly.sets['BOTTOM'].nodes:
        if round(X,4)!=abs(round(i.coordinates[0],4)) and round(Z,4)!=abs(round(i.coordinates[2],4)):

            b_BOTTOM=b_BOTTOM+[(i.coordinates[0],i.coordinates[1],i.coordinates[2],i.label)]
    b_BOTTOM.sort()
    
    nodes_number_top = len(b_TOP)
    nodes_number_bottom = len(b_BOTTOM)
    active_nodes_TB = 0
        
    rep=1
    x=0
    previous_node='not the first one' # We don't care about the first point. it will start from second one
    for i in b_TOP:
        for j in b_BOTTOM:
            if abs(i[0]-j[0])<telo*element_size and abs(i[2]-j[2])<telo*element_size  and i != previous_node:
                # Part_I and Part_II: Geometrical condition ----  Part_III: To avoid overconstrained
                active_nodes_TB +=  1
                mdb.models[modelName1].rootAssembly.Set(name='CANode-TOP-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(i[3]-1):(i[3])])
                mdb.models[modelName1].rootAssembly.Set(name='CANode-BOTTOM-'+str(rep), nodes=
                          mdb.models['Model-1'].rootAssembly.instances[instanceName1].nodes[(j[3]-1):(j[3])])
                mdb.models[modelName1].Equation(name='CConstraint-T&B-x-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 1),(1.0,'CANode-BOTTOM-'+str(rep), 1), 
                                 (0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))            
                mdb.models[modelName1].Equation(name='CConstraint-T&B-y-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 2),(1.0,'CANode-BOTTOM-'+str(rep), 2), 
                                 (0.0,'dummy1', 2),(1.0,'dummy2', 2),(0.0,'dummy3', 2)))                                    
                mdb.models[modelName1].Equation(name='CConstraint-T&B-z-'+str(x+1), 
                          terms=((-1.0,'CANode-TOP-'+str(rep), 3),(1.0,'CANode-BOTTOM-'+str(rep), 3), 
                                 (0.0,'dummy1', 3),(1.0,'dummy2', 3),(0.0,'dummy3', 3)))
                rep=rep+1
                x=x+1
                previous_node = i #We don't want to put repetitive constraint
                
                '''Note: If we delet the Y-constraint we woud have a distorsion in sides (specially in front and back).
                Try it if you wanna see!
                '''

        
    print('nodes_number_top=', nodes_number_top)
    print('nodes_number_bottom=', nodes_number_bottom)
    print('active_nodes_TB=', active_nodes_TB)
    
    if nodes_number_top < active_nodes_TB or nodes_number_bottom < active_nodes_TB:
        print('****** WARNING: OVER-CONSTRAINED ---> Top and Bottom******')

    if nodes_number_top > active_nodes_TB or nodes_number_bottom > active_nodes_TB:

        print('****** WARNING: UNDER-CONSTRAINED ---> Top and Bottom ******')
        print('There are: ', nodes_number_top - active_nodes_TB , 'nodes WITHOUT constraints')
        
    if nodes_number_top == active_nodes_TB or nodes_number_bottom == active_nodes_TB:
        print('____ Perfect CONSTRAINED ---> Top and Bottom _____')



# =============================================================================
# Constraints Shear: E12 ==> Egdes & Vertices - (to avoid overlaping)
# =============================================================================
def ConstraintsEdgesShearE12_FrRi_BaLe(last_node_FRONT_LEFT):
    modelName1='Model-1'
    ###----------------------------------------------------------------------------------     
    ## E_FrRi=E_BaLe 
    N=last_node_FRONT_LEFT         
    rep=1
    for i in range(0,N):
      	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_FrRi=E_BaLe}-'+str(i+1), 
                 terms=((-1.0,'E_FrRi-'+str(rep), 1),(1.0,'E_BaLe-'+str(rep), 1),
                        (0.0,'dummy1', 1),(0.0,'dummy2', 1),(1.0,'dummy3', 1)))
        rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_FrRi=E_BaLe}'+str(i+1), 
    	    terms=((-1.0,'E_FrRi-'+str(rep), 2),(+1.0,'E_BaLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1   
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_FrRi=E_BaLe}'+str(i+1), 
                  terms=((-1.0,'E_FrRi-'+str(rep), 3),(1.0,'E_BaLe-'+str(rep), 3),
                    (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1


def ConstraintsEdgesShearE12_BaRi_FrLe(last_node_BACK_RIGHT):       
    ###----------------------------------------------------------------------------------     
    ## E_BaRi=E_FrLe 
    modelName1='Model-1'
    N=last_node_BACK_RIGHT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaRi=E_FrLe}-'+str(i+1), 
                  terms=((-1.0,'E_BaRi-'+str(rep), 1),(1.0,'E_FrLe-'+str(rep), 1),
                    (0.0,'dummy1', 1),(0.0,'dummy2', 1),(-1.0,'dummy3', 1)))
                    
        rep=rep+1
    #-------    
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaRi=E_FrLe}'+str(i+1), 
    	    terms=((-1.0,'E_BaRi-'+str(rep), 2),(1.0,'E_FrLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1 
    #-------             
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaRi=E_FrLe}'+str(i+1), 
            terms=((-1.0,'E_BaRi-'+str(rep), 3),(1.0,'E_FrLe-'+str(rep), 3),
                   (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1 

def ConstraintsEdgesShearE12_ToLe_BoRi(last_node_TOP_LEFT):        
    # -----------------------------------------------------------------------------
    #E_ToLe=E_BoRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT    
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((1.0,'E_ToLe-'+str(rep), 1),(-1.0,'E_BoRi-'+str(rep), 1),
                (0.0,'dummy1', 1),(-1.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------  
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
    	    terms=((+1.0,'E_ToLe-'+str(rep), 2), (-1.0,'E_BoRi-'+str(rep), 2),
	       (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
        rep=rep+1
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((+1.0,'E_ToLe-'+str(rep), 3),(-1.0,'E_BoRi-'+str(rep), 3),
	        (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1


def ConstraintsEdgesShearE12_BoLe_ToRi(last_node_TOP_LEFT):           
    # -----------------------------------------------------------------------------   
    # E_BoLe=E_ToRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
		    terms=((1.0,'E_BoLe-'+str(rep), 1),(-1.0,'E_ToRi-'+str(rep), 1),
     		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------          
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
    	    terms=((1.0,'E_BoLe-'+str(rep), 2),(-1.0,'E_ToRi-'+str(rep), 2),
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1    
    #-------   
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
     	    terms=((1.0,'E_BoLe-'+str(rep), 3),(-1.0,'E_ToRi-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
    	rep=rep+1    

def ConstraintsEdgesShearE12_ToFr_BaBo(last_node_TOP_FRONT):
    ##------------------------------------------------------------------------------
    #E_ToFr=E_BaBo 
    modelName1='Model-1'       
    N=last_node_TOP_FRONT 
    rep=1
    for i in range(0,N):	
    	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	    	terms=((-1.0,'E_ToFr-'+str(rep), 1),(1.0,'E_BaBo-'+str(rep), 1),
		    (0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------   
    rep=1
    for i in range(0,N):
	    mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
    	    terms=((-1.0,'E_ToFr-'+str(rep), 2), (1.0,'E_BaBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2))); rep=rep+1
        
    #-------
    rep=1   
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	       terms=((-1.0,'E_ToFr-'+str(rep), 3),(1.0,'E_BaBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
        rep=rep+1
        

def ConstraintsEdgesShearE12_BaTo_FrBo(last_node_TOP_FRONT):
    # -----------------------------------------------------------------------------
    #E_BaTo=E_FrBo
    modelName1='Model-1'
    N=last_node_TOP_FRONT
    rep=1
    for i in range(0,N):	
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    		terms=((-1.0,'E_BaTo-'+str(rep), 1),(+1.0,'E_FrBo-'+str(rep), 1),
    		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 2),(1.0,'E_FrBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1
    rep=1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 3),(1.0,'E_FrBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
    	rep=rep+1  
    
    
def Constraints_Vertices_Shear_E12():
    modelName1='Model-1'     
    #'V_BaLeBo'-----'V_FrRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 1),(1.0, 'V_FrRiTo', 1), (-1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 2),(1.0, 'V_FrRiTo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 3),(1.0, 'V_FrRiTo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(1.0,'dummyV3', 3)))

    #'V_BaLeTo'-----'V_FrRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 1),(1.0, 'V_FrRiBo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 2),(1.0, 'V_FrRiBo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 3),(1.0, 'V_FrRiBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(1.0,'dummyV3', 3)))

    #'V_FrLeBo'----------'V_BaRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 1),(-1.0, 'V_BaRiTo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 2),(-1.0, 'V_BaRiTo' , 2), (0.0,'dummyV1', 2),(1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_BaRiTo' , 3),(-1.0, 'V_FrLeBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_FrLeTo'------'V_BaRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 1),(1.0, 'V_BaRiBo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 2),(1.0, 'V_BaRiBo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_BaRiBo', 3),(1.0, 'V_FrLeTo' , 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(1.0,'dummyV3', 3)))

   
# =============================================================================
# Constraints Shear: E13 ==> Egdes & Vertices - (to avoid overlaping)
# =============================================================================
def ConstraintsEdgesShearE13_FrRi_BaLe(last_node_FRONT_LEFT):
    modelName1='Model-1'
    ###---------------------------------------------------------------------------------     
    ## E_FrRi=E_BaLe 
    N=last_node_FRONT_LEFT         
    rep=1
    for i in range(0,N):
      	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_FrRi=E_BaLe}-'+str(i+1), 
                 terms=((-1.0,'E_FrRi-'+str(rep), 1),(1.0,'E_BaLe-'+str(rep), 1),
                        (0.0,'dummy1', 1),(0.0,'dummy2', 1),(1.0,'dummy3', 1)))
        rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_FrRi=E_BaLe}'+str(i+1), 
    	    terms=((-1.0,'E_FrRi-'+str(rep), 2),(+1.0,'E_BaLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1   
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_FrRi=E_BaLe}'+str(i+1), 
                  terms=((-1.0,'E_FrRi-'+str(rep), 3),(1.0,'E_BaLe-'+str(rep), 3),
                    (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1

def ConstraintsEdgesShearE13_BaRi_FrLe(last_node_BACK_RIGHT):       
    ###----------------------------------------------------------------------------------     
    ## E_BaRi=E_FrLe 
    modelName1='Model-1'
    N=last_node_BACK_RIGHT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaRi=E_FrLe}-'+str(i+1), 
                  terms=((-1.0,'E_BaRi-'+str(rep), 1),(1.0,'E_FrLe-'+str(rep), 1),
                    (0.0,'dummy1', 1),(0.0,'dummy2', 1),(-1.0,'dummy3', 1)))
                    
        rep=rep+1
    #-------    
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaRi=E_FrLe}'+str(i+1), 
    	    terms=((-1.0,'E_BaRi-'+str(rep), 2),(1.0,'E_FrLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1 
    #-------             
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaRi=E_FrLe}'+str(i+1), 
            terms=((-1.0,'E_BaRi-'+str(rep), 3),(1.0,'E_FrLe-'+str(rep), 3),
                   (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1 

def ConstraintsEdgesShearE13_ToLe_BoRi(last_node_TOP_LEFT):        
    # -----------------------------------------------------------------------------
    #E_ToLe=E_BoRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT    
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((1.0,'E_ToLe-'+str(rep), 1),(-1.0,'E_BoRi-'+str(rep), 1),
                (0.0,'dummy1', 1),(-1.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------  
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
    	    terms=((+1.0,'E_ToLe-'+str(rep), 2), (-1.0,'E_BoRi-'+str(rep), 2),
	       (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
        rep=rep+1
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((+1.0,'E_ToLe-'+str(rep), 3),(-1.0,'E_BoRi-'+str(rep), 3),
	        (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1

def ConstraintsEdgesShearE13_BoLe_ToRi(last_node_TOP_LEFT):           
    # -----------------------------------------------------------------------------   
    # E_BoLe=E_ToRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
		    terms=((1.0,'E_BoLe-'+str(rep), 1),(-1.0,'E_ToRi-'+str(rep), 1),
     		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------          
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
    	    terms=((1.0,'E_BoLe-'+str(rep), 2),(-1.0,'E_ToRi-'+str(rep), 2),
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1    
    #-------   
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
     	    terms=((1.0,'E_BoLe-'+str(rep), 3),(-1.0,'E_ToRi-'+str(rep), 3),
    	    (1.0,'dummy1', 3),(0.0,'dummy2', 3),(1.0,'dummy3', 3)))
    	rep=rep+1    

def ConstraintsEdgesShearE13_ToFr_BaBo(last_node_TOP_FRONT):
    ##------------------------------------------------------------------------------
    #E_ToFr=E_BaBo 
    modelName1='Model-1'       
    N=last_node_TOP_FRONT 
    rep=1
    for i in range(0,N):	
    	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	    	terms=((-1.0,'E_ToFr-'+str(rep), 1),(1.0,'E_BaBo-'+str(rep), 1),
		    (0.0,'dummy1', 1),(1.0,'dummy2', 1),(1.0,'dummy3', 1)))
        rep=rep+1
    #-------   
    rep=1
    for i in range(0,N):
	    mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
    	    terms=((-1.0,'E_ToFr-'+str(rep), 2), (1.0,'E_BaBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2))); rep=rep+1
        
    #-------
    rep=1   
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	       terms=((-1.0,'E_ToFr-'+str(rep), 3),(1.0,'E_BaBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1
        
def ConstraintsEdgesShearE13_BaTo_FrBo(last_node_TOP_FRONT):
    # -----------------------------------------------------------------------------
    #E_BaTo=E_FrBo
    modelName1='Model-1'
    N=last_node_TOP_FRONT
    rep=1
    for i in range(0,N):	
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    		terms=((-1.0,'E_BaTo-'+str(rep), 1),(+1.0,'E_FrBo-'+str(rep), 1),
    		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(-1.0,'dummy3', 1)))
    	rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 2),(1.0,'E_FrBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1
    rep=1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 3),(1.0,'E_FrBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
    	rep=rep+1  
    
    
def Constraints_Vertices_Shear_E13():
    modelName1='Model-1'     
    #'V_BaLeBo'-----'V_FrRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 1),(1.0, 'V_FrRiTo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 2),(1.0, 'V_FrRiTo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 3),(1.0, 'V_FrRiTo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_BaLeTo'-----'V_FrRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 1),(1.0, 'V_FrRiBo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 2),(1.0, 'V_FrRiBo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 3),(1.0, 'V_FrRiBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_FrLeBo'----------'V_BaRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 1),(-1.0, 'V_BaRiTo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 2),(-1.0, 'V_BaRiTo' , 2), (0.0,'dummyV1', 2),(1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_BaRiTo' , 3),(-1.0, 'V_FrLeBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_FrLeTo'------'V_BaRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 1),(1.0, 'V_BaRiBo', 1), (-1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 2),(1.0, 'V_BaRiBo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_BaRiBo', 3),(1.0, 'V_FrLeTo' , 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(1.0,'dummyV3', 3)))


    
    
# =============================================================================
# Constraints Shear: E23 ==> Egdes & Vertices - (to avoid overlaping)
# ============================================================================= 
def ConstraintsEdgesShearE23_FrRi_BaLe(last_node_FRONT_LEFT):
    modelName1='Model-1'
    ###----------------------------------------------------------------------------------     
    ## E_FrRi=E_BaLe 
    N=last_node_FRONT_LEFT         
    rep=1
    for i in range(0,N):
      	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_FrRi=E_BaLe}-'+str(i+1), 
                 terms=((-1.0,'E_FrRi-'+str(rep), 1),(1.0,'E_BaLe-'+str(rep), 1),
                        (0.0,'dummy1', 1),(0.0,'dummy2', 1),(1.0,'dummy3', 1)))
        rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_FrRi=E_BaLe}'+str(i+1), 
    	    terms=((-1.0,'E_FrRi-'+str(rep), 2),(+1.0,'E_BaLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(1.0,'dummy3', 2)))
    	rep=rep+1   
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_FrRi=E_BaLe}'+str(i+1), 
                  terms=((-1.0,'E_FrRi-'+str(rep), 3),(1.0,'E_BaLe-'+str(rep), 3),
                    (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1
        
def ConstraintsEdgesShearE23_BaRi_FrLe(last_node_BACK_RIGHT):       
    ###----------------------------------------------------------------------------------     
    ## E_BaRi=E_FrLe 
    modelName1='Model-1'
    N=last_node_BACK_RIGHT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaRi=E_FrLe}-'+str(i+1), 
                  terms=((-1.0,'E_BaRi-'+str(rep), 1),(1.0,'E_FrLe-'+str(rep), 1),
                    (0.0,'dummy1', 1),(0.0,'dummy2', 1),(-1.0,'dummy3', 1)))
                    
        rep=rep+1
    #-------    
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaRi=E_FrLe}'+str(i+1), 
    	    terms=((-1.0,'E_BaRi-'+str(rep), 2),(1.0,'E_FrLe-'+str(rep), 2), 
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(-1.0,'dummy3', 2)))
    	rep=rep+1 
    #-------             
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaRi=E_FrLe}'+str(i+1), 
            terms=((-1.0,'E_BaRi-'+str(rep), 3),(1.0,'E_FrLe-'+str(rep), 3),
                   (1.0,'dummy1', 3),(0.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1 

def ConstraintsEdgesShearE23_ToLe_BoRi(last_node_TOP_LEFT):        
    # -----------------------------------------------------------------------------
    #E_ToLe=E_BoRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT    
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((1.0,'E_ToLe-'+str(rep), 1),(-1.0,'E_BoRi-'+str(rep), 1),
                (0.0,'dummy1', 1),(-1.0,'dummy2', 1),(0.0,'dummy3', 1)))
        rep=rep+1
    #-------  
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
    	    terms=((+1.0,'E_ToLe-'+str(rep), 2), (-1.0,'E_BoRi-'+str(rep), 2),
	       (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
        rep=rep+1
    #-------      
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToLe=E_BoRi}-'+str(i+1), 
            terms=((+1.0,'E_ToLe-'+str(rep), 3),(-1.0,'E_BoRi-'+str(rep), 3),
	        (1.0,'dummy1', 3),(-1.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1

def ConstraintsEdgesShearE23_BoLe_ToRi(last_node_TOP_LEFT):           
    # -----------------------------------------------------------------------------   
    # E_BoLe=E_ToRi
    modelName1='Model-1'
    N=last_node_TOP_LEFT
    rep=1
    for i in range(0,N):
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
		    terms=((1.0,'E_BoLe-'+str(rep), 1),(-1.0,'E_ToRi-'+str(rep), 1),
     		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(0.0,'dummy3', 1)))
    	rep=rep+1
    #-------          
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
    	    terms=((1.0,'E_BoLe-'+str(rep), 2),(-1.0,'E_ToRi-'+str(rep), 2),
    	    (1.0,'dummy1', 2),(0.0,'dummy2', 2),(0.0,'dummy3', 2)))
    	rep=rep+1    
    #-------   
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BoLe=E_ToRi}-'+str(i+1), 
     	    terms=((1.0,'E_BoLe-'+str(rep), 3),(-1.0,'E_ToRi-'+str(rep), 3),
    	    (1.0,'dummy1', 3),(1.0,'dummy2', 3),(1.0,'dummy3', 3)))
    	rep=rep+1    

def ConstraintsEdgesShearE23_ToFr_BaBo(last_node_TOP_FRONT):
    ##------------------------------------------------------------------------------
    #E_ToFr=E_BaBo 
    modelName1='Model-1'       
    N=last_node_TOP_FRONT 
    rep=1
    for i in range(0,N):	
    	mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	    	terms=((-1.0,'E_ToFr-'+str(rep), 1),(1.0,'E_BaBo-'+str(rep), 1),
		    (0.0,'dummy1', 1),(1.0,'dummy2', 1),(1.0,'dummy3', 1)))
        rep=rep+1
    #-------   
    rep=1
    for i in range(0,N):
	    mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
    	    terms=((-1.0,'E_ToFr-'+str(rep), 2), (1.0,'E_BaBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(1.0,'dummy3', 2))); rep=rep+1
        
    #-------
    rep=1   
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_ToFr=E_BaBo}-'+str(i+1), 
	       terms=((-1.0,'E_ToFr-'+str(rep), 3),(1.0,'E_BaBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(1.0,'dummy2', 3),(0.0,'dummy3', 3)))
        rep=rep+1
        
def ConstraintsEdgesShearE23_BaTo_FrBo(last_node_TOP_FRONT):
    # -----------------------------------------------------------------------------
    #E_BaTo=E_FrBo
    modelName1='Model-1'
    N=last_node_TOP_FRONT
    rep=1
    for i in range(0,N):	
        mdb.models[modelName1].Equation(name='E_Constraint-x-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    		terms=((-1.0,'E_BaTo-'+str(rep), 1),(1.0,'E_FrBo-'+str(rep), 1),
    		(0.0,'dummy1', 1),(1.0,'dummy2', 1),(-1.0,'dummy3', 1)))
    	rep=rep+1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-y-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 2),(1.0,'E_FrBo-'+str(rep), 2),
    	    (0.0,'dummy1', 2),(0.0,'dummy2', 2),(-1.0,'dummy3', 2)))
    	rep=rep+1
    rep=1
    #-------
    rep=1
    for i in range(0,N):
    	mdb.models[modelName1].Equation(name='E_Constraint-z-'+'{E_BaTo=E_FrBo}-'+str(i+1), 
    	    terms=((-1.0,'E_BaTo-'+str(rep), 3),(1.0,'E_FrBo-'+str(rep), 3),
    	    (0.0,'dummy1', 3),(1.0,'dummy2', 3),(0.0,'dummy3', 3)))
    	rep=rep+1  
    
    
def Constraints_Vertices_Shear_E23():
    modelName1='Model-1'     
    #'V_BaLeBo'-----'V_FrRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 1),(1.0, 'V_FrRiTo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 2),(1.0, 'V_FrRiTo', 2), (0.0,'dummyV1', 2),(1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-1}', 
    	terms=((-1.0, 'V_BaLeBo', 3),(1.0, 'V_FrRiTo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_BaLeTo'-----'V_FrRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 1),(1.0, 'V_FrRiBo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 2),(1.0, 'V_FrRiBo', 2), (0.0,'dummyV1', 2),(1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-2}', 
    	terms=((-1.0, 'V_BaLeTo', 3),(1.0, 'V_FrRiBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(1.0,'dummyV3', 3)))

    #'V_FrLeBo'----------'V_BaRiTo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 1),(-1.0, 'V_BaRiTo', 1), (1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_FrLeBo', 2),(-1.0, 'V_BaRiTo' , 2), (0.0,'dummyV1', 2),(1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-4}', 
    	terms=((1.0, 'V_BaRiTo' , 3),(-1.0, 'V_FrLeBo', 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))

    #'V_FrLeTo'------'V_BaRiBo'
    mdb.models[modelName1].Equation(name='V_Constraint-x-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 1),(1.0, 'V_BaRiBo', 1), (-1.0,'dummyV1', 1),(0.0,'dummyV2', 1),(0.0,'dummyV3', 1)))

    mdb.models[modelName1].Equation(name='V_Constraint-y-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_FrLeTo', 2),(1.0, 'V_BaRiBo', 2), (0.0,'dummyV1', 2),(-1.0,'dummyV2', 2),(0.0,'dummyV3', 2)))

    mdb.models[modelName1].Equation(name='V_Constraint-z-'+'{T-Bo=Fr-Ba=L-R-3}', 
    	terms=((-1.0, 'V_BaRiBo', 3),(1.0, 'V_FrLeTo' , 3), (0.0,'dummyV1', 3),(0.0,'dummyV2', 3),(-1.0,'dummyV3', 3)))


        


# *****************************************************************************
' ------------------------  Boundary Conditions  -----------------------------' 
# *****************************************************************************

def BC_Axial():
    '''This functions fixed three vertices such that:
        
        'V_BaLeBo' ----> u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET
        'V_FrLeBo' ----> u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
        'V_BaRiBo' ----> u1=UNSET, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET
        
        Just for Axial E11 and E33 loading
        
  '''
    instanceName1='Part-1'

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'V_BaLeBo', region=
              mdb.models['Model-1'].rootAssembly.sets['V_BaLeBo'], u1=0.0, 
              u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)


    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'V_FrLeBo', region=
              mdb.models['Model-1'].rootAssembly.sets['V_FrLeBo'], u1=0.0, 
              u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'V_BaLeTo', region=
        mdb.models['Model-1'].rootAssembly.sets['V_BaLeTo'], u1=0.0, 
        u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'V_BaRiBo', region=
        mdb.models['Model-1'].rootAssembly.sets['V_BaRiBo'], u1=UNSET, 
        u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    
    

def BC_Shear(tole_centerP):
    
    
    '''This functions fixed the center point of the volume such that:
        
        'center_point' ----> u1=0, u2=0.0, u3=0, ur1=UNSET, ur2=UNSET, ur3=UNSET

        Note:If you have a big mesh (almost element size >  2.*c/3), this function may not work.        
  '''
    instanceName1='Part-1'
    
    Select_closet_node_to_d0(0.0,0.0,0.0,'center_point', tole_centerP)
    
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'center_point', region=
              mdb.models['Model-1'].rootAssembly.sets['center_point'], u1=0.0, 
              u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)





# *****************************************************************************
' -------------------------------  Loading  -----------------------------------'
# *****************************************************************************

def Loading_E11(E11,X,Y,Z):
    '''
    Applying E11 Strain in direction X on dummy point 1
    
    '''     
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E11-X-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy1'], 
              u1=E11*(2.*X), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)           

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E22-Y-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy2'], 
              u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E33-Z-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy3'], 
              u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)


def Loading_E22(E22,X,Y,Z):
    '''
    Applying E33 Strain in direction Z on dummy point 3
    '''
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E11-X-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy1'], 
              u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E22-Y-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy2'], 
              u1=UNSET, u2=E22*(2.*Y), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E33-Z-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy3'], 
              u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)


def Loading_E33(E33,X,Y,Z):
    '''
    Applying E33 Strain in direction Z on dummy point 3
    '''
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E11-X-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy1'], 
              u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E22-Y-Direction_FIXED', region=mdb.models['Model-1'].rootAssembly.sets['dummy2'], 
              u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E33-Z-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy3'], 
              u1=UNSET, u2=UNSET, u3=E33*(2.*Z), ur1=UNSET, ur2=UNSET, ur3=UNSET)

   
def Loading_E13(E13,X,Y,Z): 
    
    '''
    Applying E13 Strain in BOTH direction X AND Z on dummy points 1 AND 3.
    Note, in abaqus we will get gamma from odb. Here we apply E13 = 0.5 to get
    gamma_13 = 1.0, to compute one column of Cijkl from our constituve eq. with
    voight notation. 
    '''
    #----------dummy1 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E13-Z-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy1'], 
              u1=UNSET, u2=UNSET, u3=E13*(2.*X), ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummy2 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E31-X-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy3'], 
              u1=E13*(2.*Z), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    
    #----------dummyV1 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV1-U1=LZ', region=mdb.models['Model-1'].rootAssembly.sets['dummyV1'], 
              u1=-E13*(2.*Z), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummyV3 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV3-U3=LX', region=mdb.models['Model-1'].rootAssembly.sets['dummyV3'], 
              u1=UNSET, u2=UNSET, u3=E13*(2.*X), ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    
    

def Loading_E12(E12,X,Y,Z): 
    
    '''
    Applying E12 Strain in BOTH direction X AND Y on dummy points 1 AND 2.
    Note, in abaqus we will get gamma from odb. Here we apply E12 = 0.5 to get
    gamma_12 = 1.0, to compute one column of Cijkl from our constituve eq. with
    voight notation. 
    '''
    #----------dummy1 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E12-Y-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy1'], 
              u1=UNSET, u2=E12*(2.*X), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummy2 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E21-X-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy2'], 
              u1=E12*(2.*Y), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
    
    #----------dummyV1 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV1-U1=LY', region=mdb.models['Model-1'].rootAssembly.sets['dummyV1'], 
              u1=E12*(2.*Y), u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummyV2 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV2-U2=LX', region=mdb.models['Model-1'].rootAssembly.sets['dummyV2'], 
              u1=UNSET, u2=E12*(2.*X), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)


    
    

def Loading_E23(E23,X,Y,Z): 
    
    '''
    Applying E12 Strain in BOTH direction X AND Y on dummy points 1 AND 2.
    Note, in abaqus we will get gamma from odb. Here we apply E23 = 0.5 to get
    gamma_23 = 1.0, to compute one column of Cijkl from our constituve eq. with
    voight notation. 
    '''
    #----------dummy3 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E32-Y-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy3'], 
              u1=UNSET, u2=E23*(2.*Z), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummy2 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'E23-Z-Direction', region=mdb.models['Model-1'].rootAssembly.sets['dummy2'], 
              u1=UNSET, u2=UNSET, u3=E23*(2.*Y), ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummyV1 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV2-U2=LZ', region=mdb.models['Model-1'].rootAssembly.sets['dummyV2'], 
              u1=UNSET, u2=-E23*(2.*Z), u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

    #----------dummyV3 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
              distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
              'dummyV3-U3=LY', region=mdb.models['Model-1'].rootAssembly.sets['dummyV3'], 
              u1=UNSET, u2=UNSET, u3=E23*(2.*Y), ur1=UNSET, ur2=UNSET, ur3=UNSET)

    
       
# *****************************************************************************
' -------------------------------  Job  -------------------------------------' 
# *****************************************************************************   
# =============================================================================
def Job_Creat(job_name, description_note, CPU_num = 1, memory_per = 90):
    
    """
    Creat Job.
    Input:
        job_name: Name of the job
        description: job's description
        CPU_num: number of cpus using parallelization; by default CPU_num = 1
        memory_per: memory percentage usage;  by default memory_per = 90
        
    
    Hyper-parameter:
        Domains_num: number of domains corresponding to each processors.
        increas_memory: (Boolean variable) to determin whether or not auomatically increase 
        the RAM. If it's True, when software needs more memory than what you write as percentage,
        it takes more memory up to 99%. It's not recomended to use this. 
    
    Note: 
        At least 5000 DOFs are needed for each CPU core. It is not always efficient 
        to use more cores. Moreover, you need more licance (token) to use employ
        more cores. 
        
        You can increase the CPU_num up to 12
        
        It is highly recomended to consider number of domains = number of cores
        Otherwise, you would get an error. This is becasue 
            "The number of domains must be a multiple of the number of processors
            for domian level parallelization".
        
        Read Stasa page 352.
    
    
    """  
    # My defualt
    Domains_num = CPU_num
    increas_memory = True
    
    mdb.Job(name= job_name, model='Model-1', description=description_note, type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=memory_per, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=increas_memory, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT,
        numCpus=CPU_num, 
        numDomains = Domains_num,
        numGPUs=0)
    
    if CPU_num>1:
        print("Creat job with parallelization: numCpus={}=numDomains".format(CPU_num))

# =============================================================================
def Job_submit(job_name):
    'Submitting Job file'
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()
# =============================================================================

def write_job_inp(Job_name):
    mdb.jobs[Job_name].writeInput(consistencyChecking=OFF)



# _____________________________________________________________________________
' -------------------------------  F-Functions  -------------------------------' 
'                  Functions include all prevoise simulations                 ' 
# _____________________________________________________________________________

def Processing(X,Y,Z, element_size,
               telo_const_rl, telo_const_fb, telo_const_tb,
               tole_centerP,
               E11_job, E22_job, E33_job, E13_job, E12_job, E23_job):
    ''' F6_Set: Put nodes in each faces (walls), edges and vertices
    
        F7_Constraints: Put constraints between two corresponding nodes in Right/Left and Front/Back sides
        
    Parameters: 
         telo_const_(fb/rl): tolerance for putting the constraints. It may be important when 
         you have a big difference (deviation) between front/back and left/right sizes.
         
         remind ---> in my constraints function, it will be telo*element size
    '''
    
    ##Set-------------------------------------------------------------
    Functions_Set_walls(X,Y,Z)
    Functions_Set_vertices(X,Y,Z)
    (last_node_FRONT_LEFT, c_FL)=Functions_Set_Edges('FRONT', 'LEFT', 'E_FrLe',X,Y,Z)
    (last_node_FRONT_RIGHT, c_FR) = Functions_Set_Edges('RIGHT', 'FRONT', 'E_FrRi',X,Y,Z)
    (last_node_BACK_RIGHT, c_BR)=Functions_Set_Edges('RIGHT', 'BACK', 'E_BaRi',X,Y,Z)
    (last_node_BACK_LEFT, c_BL)=Functions_Set_Edges('LEFT', 'BACK', 'E_BaLe',X,Y,Z)
    
    (last_node_TOP_FRONT, c_TF)=Functions_Set_Edges('TOP', 'FRONT', 'E_ToFr',X,Y,Z)
    (last_node_TOP_BACK, c_TB)=Functions_Set_Edges('TOP', 'BACK', 'E_BaTo',X,Y,Z)
    (last_node_TOP_LEFT, c_TL)=Functions_Set_Edges('TOP', 'LEFT', 'E_ToLe',X,Y,Z)
    (last_node_TOP_RIGHT, c_TR)=Functions_Set_Edges('TOP', 'RIGHT', 'E_ToRi',X,Y,Z)
    
    (last_node_BOTTOM_FRONT, c_BoF)=Functions_Set_Edges('BOTTOM', 'FRONT', 'E_FrBo',X,Y,Z)
    (last_node_BOTTOM_BACK, c_BoB)=Functions_Set_Edges('BOTTOM', 'BACK', 'E_BaBo',X,Y,Z)
    (last_node_BOTTOM_LEFT, c_BoL)=Functions_Set_Edges('BOTTOM', 'LEFT', 'E_BoLe',X,Y,Z)
    (last_node_BOTTOM_RIGHT, c_BoR)=Functions_Set_Edges('BOTTOM', 'RIGHT', 'E_BoRi',X,Y,Z)
    
        
    Set_dummy_points(X,Y,Z, 1.2,1.2,1.2, 'dummy')
    Set_dummy_points(X,Y,Z, 1.5,1.5,1.5, 'dummyV')

    # *************************************************************************
    #  Axial: E11, E22, E33
    # *************************************************************************
    ##Constraints: Faces-------------------------------------------------------
    modelName1='Model-1'    
    Constraint_RIGHT_LEFT_Axial(element_size, telo_const_rl, X,Y,Z)   
    Constraint_FRONT_BACK_Axial(element_size, telo_const_fb, X,Y,Z)
    Constraint_TOP_BOTTOM_Axial(element_size, telo_const_tb, X,Y,Z)
    
    
    ##Constraints: Edges-------------------------------------------------------
    """
    Note: 
        For complex RVE with messy meshe, the following function may not work. 
        This is because the number of elements (or nodes) is not the same at the
        parallel edges. As a result the, for loops in it the Edge_constraint functions
        do not work correctly to put the same node in one MPC equation. 
        
        There are two ways:
            first way: Ignore
            I wrote all those funcitons in try/except structure. So it will solve 
            and compute the Cijkl regardless of the this problem. 
            
            Second way: You can use Hex-dominated Sweep/Advanced front mesh 
            algorithm for all parts (all hetergeneousities). In this way, there 
            will be (hopefully) the same number of elements in edges.  
    
    """
    
    if (last_node_FRONT_LEFT == last_node_FRONT_RIGHT and
        last_node_FRONT_LEFT ==last_node_BACK_RIGHT and
        last_node_FRONT_LEFT == last_node_BACK_LEFT and
        last_node_TOP_FRONT == last_node_TOP_BACK and
        last_node_TOP_FRONT == last_node_BOTTOM_FRONT and
        last_node_TOP_FRONT == last_node_BOTTOM_BACK and
        last_node_TOP_LEFT == last_node_TOP_RIGHT and
        last_node_TOP_LEFT == last_node_BOTTOM_LEFT and
        last_node_TOP_LEFT == last_node_BOTTOM_RIGHT):
        print('Perfect edge: Same node numbers at parallel edges')
        SameNodesAtEdges = True
        ConstraintsEdgesAxial_FrRi_BaLe(last_node_FRONT_LEFT)
        ConstraintsEdgesAxial_BaRi_FrLe(last_node_BACK_RIGHT)
        ConstraintsEdgesAxial_ToLe_BoRi(last_node_TOP_LEFT)
        ConstraintsEdgesAxial_BoLe_ToRi(last_node_TOP_LEFT)
        ConstraintsEdgesAxial_ToFr_BaBo(last_node_TOP_FRONT)
        ConstraintsEdgesAxial_BaTo_FrBo(last_node_TOP_FRONT)
        

        
    else:
        SameNodesAtEdges = False
        print('*********************************************************')
        print('Error: ')
        print('The number of nodes are not equal in the parallel edges:')
        print('------------------------------------------------------')
        print('FRONT_LEFT = {}'.format(last_node_FRONT_LEFT))
        print('FRONT_RIGHT = {}'.format(last_node_FRONT_RIGHT))
        print('BACK_RIGHT = {}'.format(last_node_BACK_RIGHT))
        print('BACK_LEFT = {}'.format(last_node_BACK_LEFT))
        print('----------------------------------------------------')
        print('TOP_FRONT = {}'.format(last_node_TOP_FRONT))
        print('TOP_BACK = {}'.format(last_node_TOP_BACK))
        print('BOTTOM_FRONT = {}'.format(last_node_BOTTOM_FRONT))
        print('BOTTOM_BACK = {}'.format(last_node_BOTTOM_BACK))
        print('----------------------------------------------------')
        print('TOP_LEFT = {}'.format(last_node_TOP_LEFT))
        print('TOP_RIGHT = {}'.format(last_node_TOP_RIGHT))
        print('BOTTOM_LEFT = {}'.format(last_node_BOTTOM_LEFT))
        print('BOTTOM_RIGHT = {}'.format(last_node_BOTTOM_RIGHT))
        print('----------------------------------------------------')
        print('Use Hex_Dominated Sweep/Advanced_Front')
        # If you want to force this code to continue, write all above contraint
        # functions wtih try/except. (But You will have distorted egdes due to 
        # lack of MPC/contraints. The results will not be accurate. Yet, they will
        # not be far from the correct values. It depends on the RVE and the mechanical
        # properties.)
        print('*********************************************************')
    
    Constraints_Vertices_Axial()
#______________________________________________________________________________
    # Note: Eij = 0.01. If you put Eij = 1.0 you will get the same results
    # I just put 0.01 to avoid large deformation since NlGeom = Off. 
    '''Job and loading for E11. 
        E11: Macro-strain in direction X. (Defult value E11=1.0)
   '''
    BC_Axial()                              # Boundary conditions for E11 
    E11=1.0 
    Loading_E11(E11,X,Y,Z)                  # Applying loading. 
    description = 'Applying E11 Strain in direction X on dummy point 1'
    Job_Creat(E11_job, description)        # Creat Job
    write_job_inp(E11_job)                  # Wrting .inp file (for linux)
    Job_submit(E11_job)                     # Submiting job
    mdb.models['Model-1'].boundaryConditions['E11-X-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E22-Y-Direction_FIXED'].suppress()
    mdb.models['Model-1'].boundaryConditions['E33-Z-Direction_FIXED'].suppress()                        
#______________________________________________________________________________
    '''Job and loading for E22. 
        E22: Macro-strain in direction Y. (Defult value E22=1.0)       
        BC: same as E11
   '''
    # E33 ------ 
    E22=1.0 
    Loading_E22(E22,X,Y,Z)                  # Applying loading.
    description = 'Applying E22 Strain in direction Y on dummy point 2'
    Job_Creat(E22_job, description)        # Creat Job
    write_job_inp(E22_job)                  # Wrting .inp file (for linux)
    Job_submit(E22_job)                     # Submiting job
    mdb.models['Model-1'].boundaryConditions['E33-Z-Direction_FIXED'].suppress()
    mdb.models['Model-1'].boundaryConditions['E22-Y-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E11-X-Direction_FIXED'].suppress()   
#                              # delete boundary conditions of this part.
#______________________________________________________________________________
    '''Job and loading for E33. 
        E33: Macro-strain in direction Z. (Defult value E33=1.0)       
        BC: same as E11
   '''
    E33=1.0   
    Loading_E33(E33,X,Y,Z)                  # Applying loading.
    description = 'Applying E33 Strain in direction Z on dummy point 3'
    Job_Creat(E33_job, description)        # Creat Job
    write_job_inp(E33_job)                  # Wrting .inp file (for linux)
    Job_submit(E33_job)                     # Submiting job
    mdb.models['Model-1'].boundaryConditions['E33-Z-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E22-Y-Direction_FIXED'].suppress()
    mdb.models['Model-1'].boundaryConditions['E11-X-Direction_FIXED'].suppress()   
#                                         #delete boundary conditions of this part.
#______________________________________________________________________________
    # ****************************************************************************
    #  Shear: E13, E12, E23
    # ****************************************************************************
    # Update Shear ------ Boundary Conditions
    mdb.models['Model-1'].boundaryConditions['V_BaLeBo'].suppress()
    mdb.models['Model-1'].boundaryConditions['V_FrLeBo'].suppress()
    mdb.models['Model-1'].boundaryConditions['V_BaLeTo'].suppress()
    mdb.models['Model-1'].boundaryConditions['V_BaRiBo'].suppress()
    BC_Shear(tole_centerP)     
    # Update Shear ------ Constraints (Faces are same for all E13, E23, E12)
    Constraint_RIGHT_LEFT_Shear(element_size, telo_const_rl, X,Y,Z)
    Constraint_FRONT_BACK_Shear(element_size, telo_const_fb, X,Y,Z)
    Constraint_TOP_BOTTOM_Shear(element_size, telo_const_tb, X,Y,Z)
    
    # Shear E13 ---------------------------------------------------------------
    if SameNodesAtEdges:
        ConstraintsEdgesShearE13_FrRi_BaLe(last_node_FRONT_LEFT)
        ConstraintsEdgesShearE13_BaRi_FrLe(last_node_BACK_RIGHT)
        ConstraintsEdgesShearE13_ToLe_BoRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE13_BoLe_ToRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE13_ToFr_BaBo(last_node_TOP_FRONT)
        ConstraintsEdgesShearE13_BaTo_FrBo(last_node_TOP_FRONT)
    
    else:    
        print('Shear E13 constraints are forced to stop.')
        
   
    '''Job and loading for E13. 
        E13: Shear Macro-strain (Defult value E13=0.5)       
        BC: We need different boundary conditions
        Constraints: We need different constraints. 
   '''  
    Constraints_Vertices_Shear_E13() 
    E13=0.5 # gamma_13 = 1.0 (in odb you will get gamma_13)
    Loading_E13(E13,X,Y,Z)                                                # Applying loading.
    description = 'Applying E13 Strain in direction X AND direction Z'
    Job_Creat(E13_job, description)                                       # Creat Job
    write_job_inp(E13_job)                                                # Wrting .inp file (for linux)
    Job_submit(E13_job)                                                   # Submiting job

    mdb.models['Model-1'].boundaryConditions['E13-Z-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E31-X-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV1-U1=LZ'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV3-U3=LX'].suppress()
    #__________________________________________________________________________
    
    # Shear E12 ---------------------------------------------------------------
    if SameNodesAtEdges:
        ConstraintsEdgesShearE12_FrRi_BaLe(last_node_FRONT_LEFT)
        ConstraintsEdgesShearE12_BaRi_FrLe(last_node_BACK_RIGHT)
        ConstraintsEdgesShearE12_ToLe_BoRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE12_BoLe_ToRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE12_ToFr_BaBo(last_node_TOP_FRONT)
        ConstraintsEdgesShearE12_BaTo_FrBo(last_node_TOP_FRONT)

    else:
        print('Shear E12 constraints are forced to stop.')
            
    '''F10_Job_Load_E12: Job and loading for E12. 
    Parametres:
        E12: Shear Macro-strain (Defult value E12=0.5)       
        BC: We need different boundary conditions
        Constraints: We need different constraints. 
        Functions --> Constraint_RIGHT_LEFT_Shear, Constraint_FRONT_BACK_Shear, Constraints_Edges_Shear
        are for updating
   '''
    Constraints_Vertices_Shear_E12() 
    E12=0.5 # gamma_12 = 1.0 (in odb you will get gamma_12)
    Loading_E12(E12,X,Y,Z)                                                # Applying loading.
    description = 'Applying E12 Strain in direction X AND direction Y'
    Job_Creat(E12_job, description)                                       # Creat Job
    write_job_inp(E12_job)                                                # Wrting .inp file (for linux)
    Job_submit(E12_job)                                                   # Submiting job

    mdb.models['Model-1'].boundaryConditions['E12-Y-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E21-X-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV1-U1=LY'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV2-U2=LX'].suppress()
    #__________________________________________________________________________
    
    # Shear E12 ---------------------------------------------------------------
#     E23 ----
    if SameNodesAtEdges:
        # Shear E23
        ConstraintsEdgesShearE23_FrRi_BaLe(last_node_FRONT_LEFT)
        ConstraintsEdgesShearE23_BaRi_FrLe(last_node_BACK_RIGHT)
        ConstraintsEdgesShearE23_ToLe_BoRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE23_BoLe_ToRi(last_node_TOP_LEFT)
        ConstraintsEdgesShearE23_ToFr_BaBo(last_node_TOP_FRONT)
        ConstraintsEdgesShearE23_BaTo_FrBo(last_node_TOP_FRONT)
    
    else:
        print('Shear E23 constraints are forced to stop.')
            
    '''F10_Job_Load_E23: Job and loading for E23. 
    Parametres:
        E23: Shear Macro-strain (Defult value E23=0.5)       
        BC: We need different boundary conditions
        Constraints: We need different constraints. 
        Functions --> Constraint_RIGHT_LEFT_Shear, Constraint_FRONT_BACK_Shear, Constraints_Edges_Shear
        are for updating
   '''
    Constraints_Vertices_Shear_E23()                                                             # Update boundary conditions
    E23=0.5 # gamma_23 = 1.0 (in odb you will get gamma_23)
    Loading_E23(E12,X,Y,Z)                                                # Applying loading.
    description = 'Applying E23 Strain in direction Y AND direction Z'
    Job_Creat(E23_job, description)        # Creat Job
    write_job_inp(E23_job)                                                # Wrting .inp file (for linux)
    Job_submit(E23_job)                                                   # Submiting job
    
    mdb.models['Model-1'].boundaryConditions['E32-Y-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['E23-Z-Direction'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV2-U2=LZ'].suppress()
    mdb.models['Model-1'].boundaryConditions['dummyV3-U3=LY'].suppress()
#______________________________________________________________________________
    

def Post_processing(E11_job, E22_job, E33_job,
                    E13_job, E12_job, E23_job, 
                    report_name):
    
    '''====================================================
    #Post-Processin ------> for Air model
    #====================================================
     REVISED 22 - 3 - 2019 BACED ON WHAT FRANCIS SAID AND PAGE 18 ELIAS PHD THESIS
     When we apply displacement we should calculate stifness
    S11     C11 C12 C13 C14 C15 C16      E11 
    S22     C21 C22 C23 C24 C25 C26      E22 
    S33  =  C31 C32 C33 C34 C36 C36 *    E33 
    S12     C41 C42 C43 C44 C46 C46    2*E12 
    S13     C51 C52 C53 C45 C55 C56    2*E13 
    S23     C61 C62 C63 C65 C65 C66    2*E23 
    
    Please read also pgae 136 of Martin's thesis. The notation in the abaqus
    is the classical Voigt notation. (Martin used Voigt notation see the appendix. 
    In writing the UMAT Subroutine, however, he had to chnage it to the classical
    notationa and then again chanage to modified notation.) This is what always bothered
    me. But I think it is correct the way that I computed the Cijkl and used the
    computed tensor in my homogenzied model. BOH (my note at 03-06-2021)'
    
    Important note:
        ABAQUS always reports shear strain as engineering shear strain, gamma_ij:
            gamma_ij = E_ij + E_ji
    '''
   

    def Stress_Strain_Ave(Job):
        'returns average stress, average strain and coresponding column of C tensor'
        Job=Job+'.odb'
        myOdb = openOdb(path=Job);
        frameRepository = myOdb.steps['Step-1'].frames;
        frameS=[]; frameE=[]; frameIVOL=[];
        frameS.insert(0,frameRepository[-1].fieldOutputs['S'].getSubset(position=INTEGRATION_POINT));
        frameE.insert(0,frameRepository[-1].fieldOutputs['E'].getSubset(position=INTEGRATION_POINT));
        frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
        Tot_Vol=0;    #Total Volume
        Tot_Stress=0; #Stress Sum
        Tot_Strain=0; #Strain Sum
        for II in range(0,len(frameS[-1].values)):
            Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
            Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;
            Tot_Strain=Tot_Strain+frameE[0].values[II].data * frameIVOL[0].values[II].data;
        Avg_Stress = Tot_Stress/Tot_Vol;
        Avg_Strain = Tot_Strain/Tot_Vol;
        C_column=[]
        for i in range(0,6):
            C_column.append(Avg_Stress[i]/Avg_Strain[i])
            
        
        return Avg_Stress, Avg_Strain 
    
    (Avg_Stress_11, Avg_Strain_11)=Stress_Strain_Ave(E11_job)
    (Avg_Stress_22, Avg_Strain_22)=Stress_Strain_Ave(E22_job)
    (Avg_Stress_33, Avg_Strain_33)=Stress_Strain_Ave(E33_job)
    
    (Avg_Stress_13, Avg_Strain_13)=Stress_Strain_Ave(E13_job)
    (Avg_Stress_12, Avg_Strain_12)=Stress_Strain_Ave(E12_job)
    (Avg_Stress_23, Avg_Strain_23)=Stress_Strain_Ave(E23_job)
    print(' ----------------------- Axial-E11 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_11))
    print('Avg_Strain={}'.format(Avg_Strain_11))
    print(' ----------------------- Axial-E22 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_22))
    print('Avg_Strain={}'.format(Avg_Strain_22))
    print(' ----------------------- Axial-E33 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_33))
    print('Avg_Strain={}'.format(Avg_Strain_33))
    print(' ----------------------- Shear-E13 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_13))
    print('Avg_Strain={}'.format(Avg_Strain_13))
    print(' ----------------------- Shear-E12 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_12))
    print('Avg_Strain={}'.format(Avg_Strain_12))
    print(' ----------------------- Shear-E23 ----------------------' )
    print('Avg_Stress={}'.format(Avg_Stress_23))
    print('Avg_Strain={}'.format(Avg_Strain_23))

    Ct=np.zeros((6,6)); #temp
    C=np.zeros((6,6));
    GPa = 1.0e9;
#    ac=3 
    ac = 9 # for convergence study 

    Ct[0,0]= (Avg_Stress_11[0]/Avg_Strain_11[0])/GPa; C[0,0]=round(Ct[0,0],ac)
    Ct[1,0]= (Avg_Stress_11[1]/Avg_Strain_11[0])/GPa; C[1,0]=round(Ct[1,0],ac)
    Ct[2,0]= (Avg_Stress_11[2]/Avg_Strain_11[0])/GPa; C[2,0]=round(Ct[2,0],ac)
    Ct[3,0]= (Avg_Stress_11[3]/Avg_Strain_11[0])/GPa; C[3,0]=round(Ct[3,0],ac)
    Ct[4,0]= (Avg_Stress_11[4]/Avg_Strain_11[0])/GPa; C[4,0]=round(Ct[4,0],ac)
    Ct[5,0]= (Avg_Stress_11[5]/Avg_Strain_11[0])/GPa; C[5,0]=round(Ct[5,0],ac)
        
    Ct[0,1]= (Avg_Stress_22[0]/Avg_Strain_22[1])/GPa; C[0,1]=round(Ct[0,1],ac)
    Ct[1,1]= (Avg_Stress_22[1]/Avg_Strain_22[1])/GPa; C[1,1]=round(Ct[1,1],ac)
    Ct[2,1]= (Avg_Stress_22[2]/Avg_Strain_22[1])/GPa; C[2,1]=round(Ct[2,1],ac) 
    Ct[3,1]= (Avg_Stress_22[3]/Avg_Strain_22[1])/GPa; C[3,1]=round(Ct[3,1],ac)
    Ct[4,1]= (Avg_Stress_22[4]/Avg_Strain_22[1])/GPa; C[4,1]=round(Ct[4,1],ac)
    Ct[5,1]= (Avg_Stress_22[5]/Avg_Strain_22[1])/GPa; C[5,1]=round(Ct[5,1],ac)
    
    Ct[0,2]= (Avg_Stress_33[0]/Avg_Strain_33[2])/GPa; C[0,2]=round(Ct[0,2],ac)
    Ct[1,2]= (Avg_Stress_33[1]/Avg_Strain_33[2])/GPa; C[1,2]=round(Ct[1,2],ac)
    Ct[2,2]= (Avg_Stress_33[2]/Avg_Strain_33[2])/GPa; C[2,2]=round(Ct[2,2],ac) 
    Ct[3,2]= (Avg_Stress_33[3]/Avg_Strain_33[2])/GPa; C[3,2]=round(Ct[3,2],ac)
    Ct[4,2]= (Avg_Stress_33[4]/Avg_Strain_33[2])/GPa; C[4,2]=round(Ct[4,2],ac)
    Ct[5,2]= (Avg_Stress_33[5]/Avg_Strain_33[2])/GPa; C[5,2]=round(Ct[5,2],ac)
      
    Ct[0,3]= (Avg_Stress_12[0]/Avg_Strain_12[3])/GPa; C[0,3]=round(Ct[0,3],ac)
    Ct[1,3]= (Avg_Stress_12[1]/Avg_Strain_12[3])/GPa; C[1,3]=round(Ct[1,3],ac)
    Ct[2,3]= (Avg_Stress_12[2]/Avg_Strain_12[3])/GPa; C[2,3]=round(Ct[2,3],ac)  
    Ct[3,3]= (Avg_Stress_12[3]/Avg_Strain_12[3])/GPa; C[3,3]=round(Ct[3,3],ac)
    Ct[4,3]= (Avg_Stress_12[4]/Avg_Strain_12[3])/GPa; C[4,3]=round(Ct[4,3],ac)
    Ct[5,3]= (Avg_Stress_12[5]/Avg_Strain_12[3])/GPa; C[5,3]=round(Ct[5,3],ac)

    Ct[0,4]= (Avg_Stress_13[0]/Avg_Strain_13[4])/GPa; C[0,4]=round(Ct[0,4],ac)
    Ct[1,4]= (Avg_Stress_13[1]/Avg_Strain_13[4])/GPa; C[1,4]=round(Ct[1,4],ac)
    Ct[2,4]= (Avg_Stress_13[2]/Avg_Strain_13[4])/GPa; C[2,4]=round(Ct[2,4],ac)
    Ct[3,4]= (Avg_Stress_13[3]/Avg_Strain_13[4])/GPa; C[3,4]=round(Ct[3,4],ac)
    Ct[4,4]= (Avg_Stress_13[4]/Avg_Strain_13[4])/GPa; C[4,4]=round(Ct[4,4],ac)
    Ct[5,4]= (Avg_Stress_13[5]/Avg_Strain_13[4])/GPa; C[5,4]=round(Ct[5,4],ac)

    Ct[0,5]= (Avg_Stress_23[0]/Avg_Strain_23[5])/GPa; C[0,5]=round(Ct[0,5],ac)
    Ct[1,5]= (Avg_Stress_23[1]/Avg_Strain_23[5])/GPa; C[1,5]=round(Ct[1,5],ac)
    Ct[2,5]= (Avg_Stress_23[2]/Avg_Strain_23[5])/GPa; C[2,5]=round(Ct[2,5],ac) 
    Ct[3,5]= (Avg_Stress_23[3]/Avg_Strain_23[5])/GPa; C[3,5]=round(Ct[3,5],ac)
    Ct[4,5]= (Avg_Stress_23[4]/Avg_Strain_23[5])/GPa; C[4,5]=round(Ct[4,5],ac)
    Ct[5,5]= (Avg_Stress_23[5]/Avg_Strain_23[5])/GPa; C[5,5]=round(Ct[5,5],ac)
    'I could write a "for", but this way is better to control stress and strain'
    

    print('''
          Constitutive equation S=C:E
    S11     C11 C12 C13 C14 C15 C16      E11 
    S22     C21 C22 C23 C24 C25 C26      E22 
    S33  =  C31 C32 C33 C34 C36 C36 *    E33 
    S12     C41 C42 C43 C44 C46 C46    2*E12 
    S13     C51 C52 C53 C45 C55 C56    2*E13 
    S23     C61 C62 C63 C65 C65 C66    2*E23 
    
    Absqus returns gamma for the shear strain (and not epsilon)
          
    ''')

    print('*\n*')
    print('------------------------------------------------------')
    print('*\n*')
    print('Stiffness tensor [GPa]')
    print(np.round(C, 3))

    np.savetxt(report_name,C,fmt='%.9f') # 9 digits for convergence study



def change_work_dir():
    os.chdir(r"D:\mcodes\homogen\abaqus_files")




# =============================================================================
# Dimentions of the cube WRT c,t,h
# =============================================================================
def f_cth_to_XZY(c,t,h):
    '''independents variables c,t,h1,h2 to---> X,Y,Z'''
    X=3.*c + np.sqrt(3)*t
    Z=np.sqrt(3)*c + t
    Y=0.5*(h)
    
    return X,Z,Y



def Test_Result():
    print('------------------------------------------------------')
    GPa=1.0e9; ac=3
    E_solid= E_void = E=1.7e9/GPa
    v_solid= v_void = v=0.4 
    print('for isotropic assumption (without void) ---> E_solid=E_void',)
    print('E_solid=',E_solid)
    print('E_void=',E_void)
    C11_iso=(E/(1+v))*(1-v)/(1-2*v)
    C12_iso=(E/(1+v))*((v)/(1-2*v))
    C44_iso=E/(2*(1+v))
    C_iso=np.zeros((6,6))
    C_iso[0,0]=C_iso[1,1]=C_iso[2,2]=round(C11_iso,ac)
    C_iso[0,1]=C_iso[0,2]=C_iso[1,0]=C_iso[1,2]=C_iso[2,0]=C_iso[2,1]=round(C12_iso,ac)
    C_iso[3,3]=C_iso[4,4]=C_iso[5,5]=round((C11_iso-C12_iso)/2.,ac)    
    print('you should have')
    print('C_iso=           [GPa]')
    print(C_iso)
    print('Ref.: Nemat-Nasser P-71')


