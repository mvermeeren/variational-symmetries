# variational-symmetries
Constructing pluri-Lagrangian stuctures using variational symmetries, as presented in

	Petrera, Vermeeren. Variational symmetries and pluri-Lagrangian structures for integrable hierarchies of PDEs, 2019.

The code consists of the following files:

	* varsym.sage
	* varsym.ipynb
	* varsym_input.sage
	* tools/2d_variational_calculus.sage
	* tools/init_continuous.sage
	* tools/output.sage

This software was developed in SageMath 8.1 by

	Mats Vermeeren
	TU Berlin
	vermeeren@math.tu-berlin.de
	
and published under the MIT License (see 'LICENCE' for more details)


RUNNING THE PROGRAM

Open the notebook varsym.ipynb in Sage (http://www.sagemath.org/). Change input as needed and load 'varsym.sage' from the notebook.

Alternatively: from the Sage command line load first 'varsym_input.sage' and then 'varsym.sage'.


OUTPUT

Output is printed in the notebook (or console).
If the parameter 'save' is set to True, the program will create log files in the 
subfolder 'logs' and save final results in the subfolder 'lagrangians'. 
If the parameter 'viewpdf' is set to True, a pdf of the output will be compiled.


ADDING MORE HIERARCHIES

Custom equations and Lagrangians can be added in the method get_input(switch) in the notebook (or in
the file 'varsym_input.sage'). A unique 'switch' should be introduced to select them.
