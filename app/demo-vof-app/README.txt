Create a gmv file with material polyhedrons:
--------------------------------

./../../build/bin/demo-vof-app3d [filename1] [filename2] [filename3]

Where [filename1] is an input binary file containing volume
fraction data for materials in mesh cells
and [filename2] is the mesh file in a format readable by Jali.

The resulting gmv output is named [filename3].

Example: 
--------
./../../build/bin/demo-vof-app3d cubic10-3M.bvf cubic10.exo cubic10-3M.gmv
or  
./../../build/bin/demo-vof-app3d cubic10-3M.bvf cubic10.exo
(this version uses [filename1].gmv as the default for [filename3])
