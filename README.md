# gji_tem_2025
Supplementary material for the GJI article

Tested with MATLAB 2022b, 2026a, Gmsh 4.15.0, and Paraview 5.13

## FEMALY

Please run the script `.../femaly/drive_paper.m` from the femaly root directory. \
The script calculates all components of the $E$ and $\dot{B}$ fields, the respective sensitivities for a single observation point and the `.pvd` file for the visualization into ParaView.

If you want to modify the script parameters, please refer to the 'Initialize.' section within `drive_paper.m`.

Run  `.../femaly/driver_plots_paper.m` to reproduce figures 5, 6, 7, and 8. \
Load `.../femaly/results/show_S.psvm` in ParaView to visualize the (pre-calculated) sensitivities in 3D and reproduce figure 4.

Default script parameters are set so that calculations can be run on a machine with 16 GB of RAM. \
To reproduce the results from the paper, change variables labeled by `% PAPER: ...`. 

MATLAB and Gmsh are required.
If the direct solver MUMPS is available, also parallel evaluation is supported.
To crosscheck the numerical solutions, empymod with Python3 must be available. 

### Setting up MUMPS
Please see https://github.com/Mathias-Scheunert/matlab-mumps-build for installation instructions.

### Setting up empymod
Please install empymod within a (mini)conda environment. \
For older MATLAB versions (2022b, not required in 2026a) set environment variables before you start MATLAB:

**bash**
```
path_gcc="<your-path>.conda/envs/<your_env>/lib/libstdc++.so.6"
export LD_PRELOAD="$path_gcc"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
```

Right after starting MATLAB, before running any scripts (or python api), run:

**MATLAB**
```
py.sys.setdlopenflags(int32(10));
```

## Convolutional Approach

First `run analytic_e_field_grid_demo.m` to produce the physical fields on a regular grid using a given loop source. \
E-field values are written into `EM_field_and_coord_empymod_5x3x2_-7_-3_srcl1.mat` where the file name contains information on the number of grid points (5x3x2), the time range -7_-3 and the source side length srcl1.

Then run `analytic_e_field_grid_demo.m` again to produce the adjoint fields on a regular grid using a small source.

Next, run `convolution_vector_3D_split_demo.m` to perform the convolution.\
It prompts for the two files created before and produces an output file named `Convolution_Results_empymod1_empymod1_grad_5x3x2_20_0_0_-7_-3.mat` where grad stands for the time derivative of the adjoint field and 20_0_0 is the location of the receiver (i.e. the location of the adjoint source). \
In this case and for simplicity, the physical and adjoint fields are generated both by empymod with a source side length of 1m. \
However, it is recommended to choose a smaller adjoint source or a dipole, particularly for fields close to the adjoint source. \
The location of the physical source is always 0_0_0. \
The program determines the overlapping spatial volume of both runs. \
The time range extends from  $10^{-7}$ to $10^{-3}$ s. \
In the demo version, the VMD option is not enabled.

The installation of empymod is a prerequisite for running `analytic_e_field_grid_demo.m`.\
Note that the E-fields can be calculated by any other routine, particularly also using FEMALY.  

The above scripts, or parts of them, were used to produce figures 3, 10, 11, 12, and 13.

## Strokkur TEM data

The Strokkur TEM data acquired on Sep 25, 2024 and presented in figures 1 and A1, can be obtained from `.../TEM_Strokkur_data/54_20240925_183542_dBdt_Z.mat`.

This file contains a structure $\dot{B}$ with fields:

| Field | Description |
| :--- | :--- |
|**UTC**| Record time stamp Universal Time|
|**Header**| Record Header including RR - repetition rate, IT - integration time, I - Tx current (device-internal default setting - not relevant), TxA - Tx coil area, RxA - Rx coil area, Gain - amplification factor, Delay - additional time shift, Equip - measuring equipment, NGates - number of gates, T0 - ramp time, Prim_feld - primary field, and Comment|
|**GateTimes**| gate times|
|**GPS_Tx**| Tx location|
|**GPS_Rx**| Rx location| 
|**Tx_Current**| Tx current (measured - relevant)|
|**X**| x-component $\dot{B}$ transient|
|**Y**| y-component $\dot{B}$ transient|
|**Z**| z-component $\dot{B}$ transient|
