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

### Setting up MUMPS
Please see https://github.com/Mathias-Scheunert/matlab-mumps-32 for install instructions.

### Setting up empymod
Please install empymod within a (mini)conda environment and set environment variables before you start Matlab:

path_gcc="<your-path>.conda/envs/<your_env>/lib/libstdc++.so.6"

export LD_PRELOAD="$path_gcc"

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

Right after starting Matlab, before running any scripts (or python api), run:

>> py.sys.setdlopenflags(int32(10));

