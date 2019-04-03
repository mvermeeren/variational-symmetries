# variational-symmetries
Constructing pluri-Lagrangian stuctures using variational symmetries, as presented in

	Petrera, Vermeeren. Variational symmetries and pluri-Lagrangian hierarchies, 2019.

The code consists of the following files:

	* varsym.sage (Main)
	* varsym_input.sage (Specify input here)
	* tools/2d_variational_calculus.sage
	* tools/init_continuous.sage
	* tools/init_discrete.sage
	* tools/output.sage

This software is written was developed in Sage 7.5.1 by

	Mats Vermeeren
	TU Berlin
	vermeeren@math.tu-berlin.de
	
and published under the MIT License (see 'LICENCE' for more details)


RUNNING THE PROGRAM
	
Execute the file 'varsym.sage' in SageMath.

All input is given in the file 'varsym_input.sage'.

The hierarchy under consideration is specified by setting 'switch'.
The dimension of the multi-time in which it is embedded is to be specified 
with 'numtimes'.
Additional options are described in the comments of the file 'varsym_input.sage'.


OUTPUT

Aside from printing in the console, the program will create log files in the 
subfolder 'logs'. If the parameter 'viewpdf' is set to True, a pdf of the
output will be compiled.

If a continuous pluri-Lagrangian structure is successfully computed, it will be 
written to a file in the subfolder 'lagrangians'


ADDING MORE HIERARCHIES

Custom equations and Lagrangians can be added in the method get_input(switch) of
the file 'varsym_input.sage'. A unique 'switch' should be introduced to select them.
