# EasyHybrid3
EasyHybrid is a free and open source graphical environment for the pDynamo3 package. It is developed using a combination of Python3, Cython3, GTK, and modern OpenGL. It utilizes the VisMol graphical engine to render 3D structures. EasyHybrid is a graphical extension of pDynamo that allows users to perform most of the basic pDynamo routines within its interface, as well as inspect, edit, and export pDynamo systems for further simulations using Python scripting in text mode.

<img width="1589" height="908" alt="image" src="https://github.com/user-attachments/assets/fd399057-fd9d-47f4-b325-92f98f270b98" />

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

# Cloning and Installation

This repository contains the latest **stable release** of EasyHybrid3.

Clone it using:

```bash
git clone --recurse-submodules https://github.com/ferbachega/EasyHybrid3.git
```

If you want the **development version**, which is updated almost weekly, you may clone the `dev` branch instead.  
Be aware that this version may contain more bugs and experimental features.

```bash
git clone --recurse-submodules -b dev https://github.com/ferbachega/EasyHybrid3.git
```

## Verify Submodules

Make sure the graphics engine was cloned correctly by checking whether the following directory exists:

```bash
../EasyHybrid3/src/graphics_engine
```

If this folder is missing, the submodules were not downloaded properly.

## Installation

Navigate to the main EasyHybrid directory:

```bash
cd ../EasyHybrid3/
```

There, you will find the installation script:

```bash
python3 ./install.py
```

> **Note:** Make sure all required dependencies are installed before running the installer.

## Running EasyHybrid

Once the graphics engine has been successfully compiled, you can launch EasyHybrid by running:

```bash
./easyhybrid.py
```


# References
pDynamo:

  The second edition of the book M. J. Field, A Practical Introduction to the Simulation of Molecular Systems published by Cambridge University Press in 2007 [Amazon UK, Amazon US, CUP].

  The paper M. J. Field, The pDynamo Library for Molecular Simulations using Hybrid Quantum Mechanical and Molecular Mechanical Potentials, J. Chem. Theo. Comp. 4, 1151–1161 (2008) [Link]. This paper provides a general introduction to pDynamo and its QC models and hybrid potentials.

  pDynamo3 paper Martin J. Field Journal of Chemical Information and Modeling 2022 62 (23), 5849-5854 DOI: 10.1021/acs.jcim.2c01239 

EasyHybrid:
  Bachega, J. F. R., Hagen, G., Sequeiros-Borja, C., Nikklas, K., Chahine, J., Timmers, L. F. M., & Field, M. J. (2025). EasyHybrid: An interactive graphical environment for quantum, classical and hybrid simulations with pDynamo3. Journal of Chemical Information and Modeling.


important links:

https://sites.google.com/view/easyhybrid/home?authuser=0  (docs and tutorials)

https://github.com/casebor/graphics_engines (engine based on OpenGL)

https://www.pdynamo.org/home (pDynamo website)

https://github.com/pdynamo/pDynamo3 

