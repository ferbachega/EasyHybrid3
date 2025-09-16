# EasyHybrid3
EasyHybrid is a free and open source graphical environment for the pDynamo3 package. It is developed using a combination of Python3, Cython3, GTK, and modern OpenGL. It utilizes the VisMol graphical engine to render 3D structures. EasyHybrid is a graphical extension of pDynamo that allows users to perform most of the basic pDynamo routines within its interface, as well as inspect, edit, and export pDynamo systems for further simulations using Python scripting in text mode.

![Screenshot from 2023-02-06 01-55-42](https://user-images.githubusercontent.com/8658227/216887855-fb1534cf-338d-4a41-8219-5cc47acca5af.png)

The interface has been designed to make it easy to determine reaction pathways in biological systems, particularly when using hybrid quantum mechanics/molecular mechanics (QC/MM) methods. EasyHybrid allows users to access a variety of pDynamo capabilities within its graphical environment, including:

Pure quantum chemistry (QC) simulations, both ab initio and semi-empirical

Pure molecular mechanics (MM) simulations using AMBER, CHARMM, or OPLS force fields

Hybrid QC/MM simulations

Single point calculations

Energy minimization

Molecular dynamics

Reaction coordinate scanning

Umbrella sampling

Reaction path calculations

# Cloning and Installing
This repository contains the latest version of EasyHybrid, including the most recent implementations and features. However, it may be more prone to occasional bugs. If you require a more stable version, consider downloading one of the stable release packages available [here](https://sites.google.com/d/1HU8-kfSypoZHc40gBQ84IEeon36i04SL/p/13CjDL5t1ceVPQXPfDzitz4GhWqYD3PBY/edit).

The installation process is straightforward, but before proceeding, please ensure that [pDynamo3](https://www.pdynamo.org/home) is already installed. Begin by cloning the repository using the following command:  
  `git clone --recurse-submodules https://github.com/ferbachega/EasyHybrid3`

Next, navigate to the ../EasyHybrid3/src/graphics_engine folder. Here, you'll find a script for compiling the graphics engine. Run the install.sh script from within this folder:
  `cd src/graphics_engine`
and 
  `./install.sh`

If the graphics_engine folder appears to be empty, you will need to clone it. To do so, execute the following command from within the src folder:
  `git clone https://github.com/casebor/graphics_engine` 
  
Once you have successfully compiled the graphics engine, simply run the easyhybrid.py file located in the EasyHybrid base directory:  
  `./easyhybrid.py` 

For more detailed information, please visit the [EasyHybrid documentation](https://sites.google.com/view/easyhybrid/user-guide?authuser=1).


# References
pDynamo:

  The second edition of the book M. J. Field, A Practical Introduction to the Simulation of Molecular Systems published by Cambridge University Press in 2007 [Amazon UK, Amazon US, CUP].

  The paper M. J. Field, The pDynamo Library for Molecular Simulations using Hybrid Quantum Mechanical and Molecular Mechanical Potentials, J. Chem. Theo. Comp. 4, 1151–1161 (2008) [Link]. This paper provides a general introduction to pDynamo and its QC models and hybrid potentials.

  pDynamo3 paper Martin J. Field Journal of Chemical Information and Modeling 2022 62 (23), 5849-5854 DOI: 10.1021/acs.jcim.2c01239 

EasyHybrid:
  The paper J. F. R. Bachega, L. F. S. M. Timmers, L. Assirati, L. B. Bachega, M. J. Field, T. Wymore. J. Comput. Chem. 2013, 34, 2190-2196. DOI: 10.1002/jcc.23346



important links:
https://github.com/casebor/graphics_engines (engine based on OpenGL)

https://www.pdynamo.org/home (pDynamo website)

https://github.com/pdynamo/pDynamo3 

