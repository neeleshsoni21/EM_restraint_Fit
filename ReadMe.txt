#To run the following test case:

#-------------------------------------
# Working code
#-------------------------------------


cd EM_Restraint_Bug
python src/TestEM_restraint.py

The code should produce the rmf trajectory.


#-------------------------------------
# Bug Reproduce
#-------------------------------------
Change line 41 in TestEM_restraint.py:

FROM:
	
	density_residues_per_component=10,

TO:

	density_residues_per_component=1,


The code then produces following exception:
	
	ValueError: Expected n_samples >= n_components but got n_components = 745, n_samples = 744

744 is the number of residues in Nup85.
