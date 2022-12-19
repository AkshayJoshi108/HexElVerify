Purpose of folder:
To validate the deformation gradients used by Abaqus. 

Method adopted:
An abaqus simulation is run for a geometry with both- single element and multiple elements (hex elements C38DR).
By including the following line at the end of the input file (before *end) we can get abaqus to output the step-wise deformation gradients at 
every integration point to a '%Jobname'.dat file.

These deformation gradients are then compared with what we can infer from the nodal displacement data using getNodalDefGrads.py
The above python code outputs the infered deformation gradients to a '%jobname'_DefGrads.csv file. 
The python script is run using abaqus's run-script command. 

Results:
The results for defgrads from Abaqus and the displacement-infered defgrads match well with relative reside of 1e-4.

