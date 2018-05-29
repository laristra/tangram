Create a gmv file with material polyhedrons:
--------------------------------

./../../build/bin/demo-vof-app3d [nx] [ny] [nz] [filename1] [filename2]

Where [filename1] is an input binary file containing volume
fraction data for materials in mesh cells.

The resulting gmv output is named [filename2]
[nx], [ny], and [nz] are the number of cells in x, y, and z directions

Example: 
--------

./../../build/bin/demo-vof-app3d 20 20 10 eye_20x20x10_vfracs.txt eye_20x20x10.gmv
or  
./../../build/bin/demo-vof-app3d 20 20 10 eye_20x20x10_vfracs.txt
(this version uses [filename1].gmv as the default for [filename2])