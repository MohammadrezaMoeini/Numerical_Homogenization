# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:02:54 2022

@author: momoe
"""

import os

# =============================================================================
# Run the required scripts
# =============================================================================
execfile('D:/mpcodes/Numerical_Homogenization/Compute_Cij.py', __main__.__dict__)
execfile('D:/mpcodes/Numerical_Homogenization/RVE_Honeycomb.py', __main__.__dict__)


# =============================================================================
# The main file for the Example 01
# =============================================================================
def Run_Example01():
    
    # ---------------------------------------------------------------------
    # Create and Mesh the RVE 
    # ---------------------------------------------------------------------
    # Geometrical parameters 
    c=2.0; t=1.0; h=4.0 # for the relative density of rho=40%
    (X,Z,Y) = f_cth_to_XZY(c,t,h)

    # Mechanical properties 
    E_Solid = 2.5e9; v_Solid=0.4; name_Solid = 'PLA' 
    E_Air =  2.5; v_Air=0.4; name_Air = 'Void'  

    # Meshing parameters 
    element_size=h/4.
    dF=0.1
    min_size_fac=1.0
    edge_fb=element_size
    edge_rl=element_size 
    
    # Geometry
    F1_Geometry(c,t,h)
    F2_Materials(E_Solid, v_Solid, name_Solid, E_Air, v_Air, name_Air)
    F3_Assembly(X,Y,Z)

    # Step: static 
    Step()

    # Mesh
    edge_seed_fb(edge_fb, dF)
    edge_seed_rl(edge_rl, dF)
    partition_cell()
    Mesh_HEX_DOMINATED_Q_Advance_front(element_size,min_size_fac, dF)
       
    # ---------------------------------------------------------------------
    # Homogenization-------------------------------------------------------
    # ---------------------------------------------------------------------
    # Constraints to find the corresponding nodes 
    telo_const_rl = 0.25 
    telo_const_fb = 0.1    
    telo_const_tb = 0.1              
    telo_const_rl_shear = 0.25
    telo_const_fb_shear = 0.1
    tole_centerP = 0.5 
    
    # Jobs 
    E11_job = 'E11_Sim'  
    E22_job = 'E22_Sim'
    E33_job = 'E33_Sim'
    E13_job = 'E13_Sim'  
    E12_job = 'E12_Sim'
    E23_job = 'E23_Sim'
    report_name = 'Cij_PLA' +'_Ex01Honeycomb' + '.txt' 

    # Run the simulaitons
    Processing(X,Y,Z, element_size,
               telo_const_rl, telo_const_fb, telo_const_tb,
               tole_centerP,
               E11_job, E22_job, E33_job,
               E13_job, E12_job, E23_job)

    # Extract the stress and strain and compute the Cijkl
    Post_processing(E11_job, E22_job, E33_job,
                    E13_job, E12_job, E23_job,
                    report_name)
    


# =============================================================================
# Change the directory 
# =============================================================================
directory = r"D:\mpcodes\Numerical_Homogenization\Examples\Example01\AabaqusFiles"
change_work_dir(directory)

Run_Example01()