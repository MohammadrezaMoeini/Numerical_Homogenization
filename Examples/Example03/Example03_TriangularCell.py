# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:17:04 2022

@author: momoe
"""

import numpy as np
import os

# =============================================================================
# Run the required scripts
# =============================================================================
execfile('D:/mpcodes/Numerical_Homogenization/Compute_Cij.py', __main__.__dict__)
execfile('D:/mpcodes/Numerical_Homogenization/RVE_Triangle.py', __main__.__dict__)

# =============================================================================
# Simualtion function  
# =============================================================================  
def Compute_props_Triangle(a, t, element_size, tol_centP_list,
                           SimName, report_name): 
    print('*****************************************************************')
    print("            Simulation {}".format(SimName))
    print('*****************************************************************')
    # Start simulation
    Mdb()  
    # Change directory --------------------------------------------------------
    mainfolder = r"D:\mpcodes\Numerical_Homogenization\Examples\Example03\AbaqusFiles"
    os.chdir(mainfolder) # It is better to change dir. inside the for loop
    # Material properties------------------------------------------------------
    # Solid
    E_Solid = 2.05563e9; v_Solid=0.37; 
    name_Solid = 'Nylon12'; sections_Solid = 'solid_section'
    # Void
    E_Void =  2.05563; v_Void=0.37;
    name_Void = 'Void'; sections_Void = 'void_section'
    # Mesh parameters ---------------------------------------------------------
    min_size_fac = 0.1; dF = 0.1
    # =========================================================================
    # Run Simulation
    # =========================================================================
    p1, p2, p3, p4, p5, p6, p7, p8, p9, q1, q2, q3, q4, q5, q6, q7, q8, q9, X, Z, pp5, pp9, qq5, qq9 = compute_points(a, t)
    Y = Z
    Cube(X,Y,Z)
    Cut_traingles(X, Y, Z,
               p1, p2, p3, p4, p5, p6, p7, p8, p9,
               q1, q2, q3, q4, q5, q6, q7, q8, q9)
    Fill_traingles(X, Y, Z,
               p1, p2, p3, p4, p5, p6, p7, p8, p9,
               q1, q2, q3, q4, q5, q6, q7, q8, q9)
    Material_Section_Solid(E_Solid, v_Solid, name_Solid, sections_Solid)
    Material_Section_Void(E_Void, v_Void, name_Void, sections_Void)
    Material_color()
    Assembly()                 # Assembly the model
    Move_Ref_to(0.0, 0.0, -Z)  # Gives the origin of coordinates at the center of model
    Step()
    Mesh(element_size,min_size_fac, dF)
    # ---------------------------------------------------------------------
    # Homogenization-------------------------------------------------------
    # ---------------------------------------------------------------------
    #_______________Constrain'
    telo_const_rl = 0.25 #0.25 #for air 0.1 is good -----> right and lef
    telo_const_fb = 0.1  #0.1 #for air 0.1 is good -----> front and back  
    telo_const_tb = 0.1              
    tole_centerP = tol_centP_list[i]            
    #_______________Jobs'
    #_______________Jobs'
    E11_job = 'E11_Sim_' + SimName
    E22_job = 'E22_Sim_' + SimName
    E33_job = 'E33_Sim_' + SimName
    E13_job = 'E13_Sim_' + SimName
    E12_job = 'E12_Sim_' + SimName
    E23_job = 'E23_Sim_' + SimName
     
    Processing(X,Y,Z, element_size,
               telo_const_rl, telo_const_fb, telo_const_tb,
               tole_centerP,
               E11_job, E22_job, E33_job,
               E13_job, E12_job, E23_job)
    
    Post_processing(E11_job, E22_job, E33_job,
                    E13_job, E12_job, E23_job,
                    report_name)



# =============================================================================
# Run function 
# =============================================================================
t_list = np.array([1.6801]) # for a=10.0
# correspond to 
rho_list = np.array([40.0])
element_list = np.array([1.0])
tol_centP_list = np.array([2.0])   
        
for i in range(0, len(rho_list)):
    a = 10.0
    t = t_list[i]  
    element_size = element_list[i];
    SimName = 'TriangleCell_rho' + str(int(rho_list[i]))  
    report_name = 'Cij_Nylon12_Tr_rho'+str(int(rho_list[i]))+'_R05072022.txt'
    Compute_props_Triangle(a, t, element_size, tol_centP_list, 
                         SimName, report_name)
    
    
    

