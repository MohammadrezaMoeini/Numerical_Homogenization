# Introduction
This repository includes python scripts using in Abaqus to compute the effective
properties of a heterogeneous materials with linearly elastic behaviour. 

# Background
## Numerical Homogenization
The homogenization technique provides the effective properties of a heterogeneous material 
based on the mechanical properties of its micro-structure. In numerical homogenization,
the effective strain tensor is computed by applying six orthogonal and 
known macro-strain tensors ($E_{ij}$) over the Representative Volume Element 
(RVE) of the material:
$$\mathbf{\Sigma} = \mathbf{\tilde{C}} \mathbf{E}$$

The macro-stress ($\Sigma_{ij}$) is computed by a volume integration over the
whole RVE: $ \mathbf{\Sigma}=\frac{1}{V} \int \boldsymbol{\sigma}(\mathbf{x})\text{d}V$.
Then, each column of the sttiffness tenosr $\tilde{C}_{ijkl}$ is computed by each 
loading case (i.e. those macro-strain)


## Periodic Boundary conditions
Periodic Boundary Conditions (PBC) should be applied over the parallel faces of the cubic
RVE. In this way, between each pair of corresponding nodes on the opposite faces
(the nodes sharing the same in-plane coordinates for two surfaces having the same normal)
the displacement PBC is satisfied when:
$$\mathbf{u}(\mathbf{x_2})-\mathbf{u}(\mathbf{x_1})=\mathbf{E}\cdot(\mathbf{x_2}-\mathbf{x_1})$$

where x1 and x2 are location vectors. Moreover,
$\mathbf{u}(\mathbf{x_1})$ and $\mathbf{u}(\mathbf{x_2})$ are the displacement
vectors of two corresponding nodes. 



# Quick Start
This is a very simple example for the honeycomb RVE.
Open Abaqus and run the following codes: File/Run Script...

## Example 01: Honeycomb cell 
Run the Example01_Honeycomb.py (don't forget to change the directory into your disk)

The script will;
1. Create one honeycomb cell as the RVE.
2. Mesh the RVE for the given element size
3. Run 6 simulations applying the 6 orthogonal macro-strian tensors
4. Extract the micro-stress (and also volume elements and micro-strain) from odb files
5. Compute the average strain (i.e. the macro-stress)
6. Compute the effective properties and save the Cijkl in the directory

**Notes**: 
- The verification and validation of the computed effective properties of honeycomb
cells are published in:
Moeini, Mohammadreza, Mickael Begon, and Martin Lévesque. 
"Numerical homogenization of a linearly elastic honeycomb lattice structure
and comparison with analytical and experimental results." Mechanics of Materials 167 (2022): 104210.


## Example 02: Square cell 
Run the Example02_SquareCell.py 
This script is the same as Example #1 but of the square cell. It also consists
of a for loop to run for different geometrical parameters or element size (for
the convergence study).

## Example 03: Triangular cell 
Run the Example03_TriangularCell.py 
This script is also the same as the previuous lattice cells. It is for triangular
cell having linearly elastic properties. 

**Notes**: 
- Comparing the effective properties of the square and triangular cells
(Example #2 and #3) are published in: 
Moeini, Mohammadreza, Anne-Laure Ménard, Lingyu Yue, Maryam Hajizadeh, 
Mickael Begon, and Martin Lévesque. "Computationally efficient model to 
predict the deformations of a cellular foot orthotic." Computers in Biology and 
Medicine (2022): 105532. 


# Possible issues 
## Abaqus version
The cods were wrriten for Abaqus 6.14, and might have some issues for the higher version. 

## Meshing 
If the RVE is geometrically complex (RVE with random fiber distribution), it is 
possible that the mesh is not identical in the opposite faces. As a result of that, the 
algorithm might no be able to find the corresponding node. Or it would not find 
the correct matches. 

## The center point
A node close to the centerpoint of the RVE was defined and fixed to avoid the 
rigid body motion. In the coarse meshes (big element sizes), there might be an error
to define this center node. In this case, you can either reduce the element size
or increase the tolerence `tole_centerP`. 



# Dependencies
- Numpy 
- Abaqus 6.14


# Citation
- @article{moeini2022numerical,
  title={Numerical homogenization of a linearly elastic honeycomb lattice structure and comparison with analytical and experimental results},
  author={Moeini, Mohammadreza and Begon, Mickael and L{\'e}vesque, Martin},
  journal={Mechanics of Materials},
  volume={167},
  pages={104210},
  year={2022},
  publisher={Elsevier}
}

For the homogenization of the square and triangular cells, please cite:
- @article{moeini2022computationally,
  title={Computationally efficient model to predict the deformations of a cellular foot orthotic},
  author={Moeini, Mohammadreza and M{\'e}nard, Anne-Laure and Yue, Lingyu and Hajizadeh, Maryam and Begon, Mickael and L{\'e}vesque, Martin},
  journal={Computers in Biology and Medicine},
  pages={105532},
  year={2022},
  publisher={Elsevier}
}

 
