To test run command:
../xyz2cub.py test.xyz

Warning:
In testing this code .xyz files have been found that do not increment xyz in the 
correct order. If the grid dimensions of the .cub file are incorrect check if this is 
the case. If not, the next most likely points of failure are lack of a new line
ending the .xyz file and the rel\_tol variable for the isClose calls being a bad size. 
