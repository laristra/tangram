Create a gmv file with material polyhedrons:
--------------------------------

./../../build/bin/demo-vof-app3d [nx] [ny] [nz] [filename1] [filename2]

Where [filename1] is an input text file containing volume
fraction (vf) data for each material (mat) in each cell, and is 
formatted as follows:

[number of mats in problem]
[vf for cell 1, mat 1] [vf for cell 1, mat 2] ...
[vf for cell 2, mat 1] [vf for cell 2, mat 2] ...
.
.
.

The resulting gmv output is named [filename2]
[nx], [ny], and [nz] are the number of cells in x, y, and z directions

Example: 
--------

./../../build/bin/demo-vof-app3d 20 20 10 eye_20x20x10_vfracs.txt eye_20x20x10.gmv
or  
./../../build/bin/demo-vof-app3d 20 20 10 eye_20x20x10_vfracs.txt
(this version uses [filename1].gmv as the default for [filename2])
or
./../../build/bin/demo-vof-app3d
(this version uses nx=ny=nz=6, 
3d_diamond_6x6x6_vfracs.txt as the default for [filename1]
and 3d_diamond_6x6x6.gmv as the default for [filename2])
