# TurbomoleScripts
## OutputParsing
All scripts that can be used for analyzing the output of a Turbomole 
calculation from a release version of Turbomole. 
Run [scriptname] -h on any script for description and usage information. 

Summary of scripts:
- backUp:      A Bash script that backs up files from local and remote 
               hosts. Makes an extra timestamped version of a backed up
               file any time the file changes. 
- getRotConst: Returns the rotational constants for the molecule 
               specified within a Turbomole format coord file. 
- xyz2cub:     Converts a Turbomole format xyz file to a cube file. A
               Turbomole format xyz file is NOT a coordinate file. It
               holds property values on a grid. 
