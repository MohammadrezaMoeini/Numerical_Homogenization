# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 10:31:26 2022

@author: momoe

Example #2: Hmogenization of square cell having linearly elastic material properties.
In this example, a for loop was added to the code to run the simulation for several input
parameters such different geometrical parameters or element size (for convergence
study).

"""
import numpy as np
import os

# =============================================================================
# Run the required scripts
# =============================================================================
execfile('D:/mpcodes/Numerical_Homogenization/Compute_Cij.py', __main__.__dict__)
execfile('D:/mpcodes/Numerical_Homogenization/RVE_Square.py', __main__.__dict__)


# =============================================================================
# Main Simulation function  
# =============================================================================  
def Run_Example02(lb, t, element_size, tol_centP_list,
                         SimName, report_name): 
    print('*****************************************************************')
    print("                 Simulation {}".format(SimName))
    print('*****************************************************************')
    # Start simulation
    Mdb()  
    # Change directory --------------------------------------------------------  
    mainfolder = r"D:\mpcodes\Numerical_Homogenization\Examples\Example02\AbaqusFiles"
    os.chdir(mainfolder) # It's better to be inside the for loop
    # Material properties------------------------------------------------------
    # Solid
    E_Solid = 2.05563e9; v_Solid=0.37; 
    name_Solid = 'Nylon12'; sections_Solid = 'solid_section'
    # Void
    E_Void =  2.05563; v_Void=0.37;
    name_Void = 'Void'; sections_Void = 'void_section'
    # Mesh parameters ---------------------------------------------------------
    min_size_fac = 0.1; dF = 0.1
    # Compute the required points to create the cube---------------------------
    a = lb + t
    X = a/2.
    Y = a/2.
    Z = a/2.   
    p1, p2, p3, p4 = compute_points(X, Y, Z, t)
    # =========================================================================
    # Run Simulation
    # =========================================================================
    Cube(X,Y,Z)
    Cut_square(X, Y, Z, p1, p2, p3, p4)
    Fill_square(X, Y, Z, t, p1, p2, p3, p4)
    Material_Section_Solid(E_Solid, v_Solid, name_Solid, sections_Solid)
    Material_Section_Void(E_Void, v_Void, name_Void, sections_Void)
    Material_color()
    Assembly()                 # Assembly the model
    Move_Ref_to(0.0, 0.0, -Z)  # Gives the origin of coordinates at the center
    Step()
    Mesh(element_size,min_size_fac, dF)
    make_it_Quadratic()
    # ---------------------------------------------------------------------
    # Homogenization-------------------------------------------------------
    # ---------------------------------------------------------------------
    #_______________Constrain'
    telo_const_rl = 0.01 # previous value 0.25 
    telo_const_fb = 0.01 # previous value 0.1 
    telo_const_tb = 0.01  # previous value 0.1             
    tole_centerP = tol_centP_list[i] # default = 0.5             
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
# Geometrical Parameters -------- TEST
lb_list = np.array([6.0])
t_list = np.array([1.746])
rho_list = np.array([40.])
element_list_FineMesh = np.array([1.5])

tol_centP_list = np.array([2.0]) 

for i in range(0, len(rho_list)):
    lb = lb_list[i]
    t = t_list[i]  
    element_size = element_list_FineMesh[i];

    SimName = 'Sq_CoMesh10_rho' + str(int(rho_list[i]))  + 'R05072022' 
    report_name = 'Cij_N12_Sq_CoMesh10_rho'+str(int(rho_list[i]))+'R05072022.txt'
    Run_Example02(lb, t, element_size, tol_centP_list, 
                         SimName, report_name)
    
    
    
    
    
    