This demo gives an example of how the convolutional sensitivities are calculated.

First run analytic_e_field_grid_demo.m to produce the physical fields on a regular grid using a given loop source. 
E-field values are written into EM_field_and_coord_empymod_5x3x2_-7_-3_srcl1.mat where the 
file name contains information on the number of grid points (5x3x2), the time range -7_-3 and the source side length srcl1.

Then run analytic_e_field_grid_demo.m again to produce the adjoint fields on a regular grid using a small source or a VMD.

Next, run convolution_vector_3D_split_demo.m to perform the convolution. It prompts for the two files created before and produces 
an output file named Convolution_Results_empymod1_empymod1_grad_5x3x2_20_0_0_-7_-3.mat where grad stands for the time derivative of the 
adjoint field and 20_0_0 is the location of receiver (i.e. the location of the adjoint source). In this case and for simplicity, the physical and adjoint fields 
are generated both by empymod with a source side length of 1m. However, it is recommended to choose a smaller adjoint source or a dipole. The location of the physical source is always 0_0_0. The program determines the overlapping spatial volume of both runs. The time range extends from 1e-7 to 1e-3 s.

The installation of empymod is a prerequisite for running analytic_e_field_grid_demo.m. Note that the E-fields can be calculated by any other routine, particularly also using FEMALY. 