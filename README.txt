This is a basic finite element (FE) script for 3D analysis of structural components, called here **FE Project**. The code was written with the same architecture of Abaqus, so that the topology optimization script BESO (Bi-directional Structural Topology Optimization) by Zuo and Xie also works here

# Example model

For an example of how to create an input file (.rcae), submit, save results and plot results, check the file "**Example_model_creation.py**". This is a classical MBB beam problem where displacement, Von Mises stress, and compliance field are calculated.

The code makes use of functions that are highly related to Abaqus scripting functions to create the model. At the moment, only a hexagonal linear mesh is available. The mesh can be manually created (elements connectivity and nodes coordinates), but in the example, it was created with a simplistic function for cubic structures.

Plotting functions of field variables are also available through the use of **matplotlib**.

# BESO

The original script for Abaqus from BESO was here slightly adapted to work with FE Project. Basically, the packages to be imported were changed and the output of the results. In this version, the output files are saved in a folder called "Output files", as well as a picture of the current status of the evolution. At the end of the process, a .txt file is generated with the history results. For more information about this script for BESO, check out: Zuo ZH, Xie YM (2015) A simple and compact Python code for complex 3D topology optimization.