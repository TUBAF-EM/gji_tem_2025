# gji_tem_2025
Supplementary material for the GJI article

Tested with Matlab 2022b and Gmsh 4.15.0

## FEMALY

Please run the script .../femaly/+app_tem/+examples/drive_paper.m from the femaly root directory.
The script calculates all components of the E and dB/dt field as well as the respective sensitivities for a single observation point.

Matlab 2022b and Gmsh are required.
If the direct solver MUMPS is available, also parallel evaluation is supported.
To crosscheck the numerical solutions, empymod with Python3 needs to be available. 

If you want to modify the script parameters, please refer to the 'Initialize.' section within drive_paper.m.

You can load the stage file .../femaly/results/show_S.psvm into ParaView to visualize the pre-calculated sensitivities in 3D.

### Setting up MUMPS
Please see https://github.com/Mathias-Scheunert/matlab-mumps-32 for install instructions.

### Setting up empymod
Please install empymod within a (mini)conda environment and set environment variables before you start Matlab:

path_gcc="<your-path>.conda/envs/<your_env>/lib/libstdc++.so.6"

export LD_PRELOAD="$path_gcc"

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

Right after starting Matlab, before running any scripts (or python api), run:

>> py.sys.setdlopenflags(int32(10));

## Convolutional Approach

This demo gives an example of how the convolutional sensitivities are calculated.

First run analytic_e_field_grid.m to produce the physical fields on a regular grid using a given loop source. 
E-field values are written into EM_field_and_coord_empymod_5x3x2_-7_-3_srcl1.mat where the 
file name contains information on the number of grid points (5x3x2), the time range -7_-3 and the source side length srcl1.

Then run analytic_e_field_grid_demo.m again to produce the adjoint fields on a regular grid using a small source or a VMD.

Next, run convolution_vector_3D_split_demo.m to perform the convolution. It prompts for the two files created before and produces 
an output file named Convolution_Results_empymod1_empymod1_grad_5x3x2_20_0_0_-7_-3.mat where grad stands for the time derivative of the 
adjoint field and 20_0_0 is the location of receiver (i.e. the location of the adjoint source). In this case and for simplicity, the physical and adjoint fields 
are generated both by empymod with a source side length of 1m. However, it is recommended to choose a smaller adjoint source or a dipole. The location of the physical source is always 0_0_0. The program determines the overlapping spatial volume of both runs. The time range extends from 1e-7 to 1e-3 s.

The installation of empymod is a prerequisite for running analytic_e_field_grid_demo.m. Note that the E-fields can be calculated by any other routine, particularly also using FEMALY. 
