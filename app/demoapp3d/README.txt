Create data for python plotting:
--------------------------------

./../../build/bin/demoapp3d [filename1] > [filename2]

Plot data:
----------

python plot.py [filename2]

Where [filename1] is an input text file containing volume
fraction (vf) data for each material (mat) in each cell, and is 
formatted as follows:

[number of mats in problem]
[vf for cell 1, mat 1], [vf for cell 1, mat 2], ...
[vf for cell 2, mat 1], [vf for cell 2, mat 2], ...
.
.
.

The resulting python plot is named [filename2].png

Example: 
--------

/../../build/bin/demoapp3d 3d_diamond_6x6x6_vfracs.txt > 3ddiamondout.txt
or  
/../../build/bin/demoapp3d > 3ddiamondout.txt
(this version uses 3d_diamond_6x6x6_vfracs.txt as the default for [filename1])

python plot.py 3ddiamondout.txt

The resulting python plot is named 3ddiamondout.txt.png
