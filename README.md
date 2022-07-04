# Numerical Homogenization
Homogenization technique provides the effective properties effective properties 
of a hetergenous material based on the mechanical properties of its microstrcut.
In numerical homogenization, by applying six orthogonal strain tensors over the 
*Representetive Volume Element* (RVE) of the material, the sttiffness tenosr Cijkl
is computed by a vlume integration of the micro-stress and micro-strain. 
The microstrures should consist of all hetergenousiry (crack, fiber, prosity, etc)


Numerical homogenization provides the effective properties of a hetergenous
material. $\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$

$$\sigma_{ij}=C_{ijkl}\varepsilon_{kl}$$



# Quick Start
This is very simple example for the honeycomb RVE. 

# Lattice structures
Example for the honeycmb cell, square and triangular cells



# Dependencies
Numpy 
Abaqus 6.14


# Citation
Honeycomb homogenization:
@article{moeini2022numerical,
  title={Numerical homogenization of a linearly elastic honeycomb lattice structure and comparison with analytical and experimental results},
  author={Moeini, Mohammadreza and Begon, Mickael and L{\'e}vesque, Martin},
  journal={Mechanics of Materials},
  volume={167},
  pages={104210},
  year={2022},
  publisher={Elsevier}
}

For the square and triangular homoegenization 
@article{moeini2022computationally,
  title={Computationally efficient model to predict the deformations of a cellular foot orthotic},
  author={Moeini, Mohammadreza and M{\'e}nard, Anne-Laure and Yue, Lingyu and Hajizadeh, Maryam and Begon, Mickael and L{\'e}vesque, Martin},
  journal={Computers in Biology and Medicine},
  pages={105532},
  year={2022},
  publisher={Elsevier}
}
