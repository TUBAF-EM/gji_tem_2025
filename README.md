# gji_tem_2025
Supplementary material for the GJI article

## FEMALY
Please run .../femaly/+app_tem/+examples/drive_paper.m from the .../femaly main folder.
The script calculates all three components of the electrical and magnetic field and the respective sensitivities for a single observation point.
The model geometry corresponds to that, described in the paper.

Please note, that time discretization via
```
t_steps_decade = 42;
```
and / or spatial discretization via
``` 
size_at_wr = src_length / 2;
size_at_pt = 5e-1;
size_at_box = size_at_pt * 3;
ref_global = 0;
```
need to be adjusted in order to ensure appropriate low discretization errors in space and time.

### MUMPS
Please see https://github.com/Mathias-Scheunert/matlab-mumps-32 for install instructions.

