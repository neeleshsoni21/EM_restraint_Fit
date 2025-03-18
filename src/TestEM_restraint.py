

import IMP
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.pmi
import IMP.pmi.dof
import IMP.pmi.macros
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Parse commandline Inputs
parser_obj = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
parser_obj.add_argument("-o", "--output_dir", default="output", help="Output files directory")
parser_obj.add_argument("-d", "--data_dir", default="data", help="Data files directory")

args = parser_obj.parse_args()
data_path = os.path.abspath(args.data_dir)
out_path = os.path.abspath(args.output_dir)

include_gaussian_em_restraint=True

gmmfilename = os.path.join(data_path,'gmms/nup85_Simulated_Map_ng100.txt')
nup85_pdb_fn=os.path.join(data_path,'pdbs/nup85_AF.pdb')
nup85_fastafile = os.path.join(data_path,'fasta/nup85.fasta')

nup85_seqs = IMP.pmi.topology.Sequences(nup85_fastafile)

bskt_mols =[]

mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

nup85_mol0 = st.create_molecule('nup85', sequence=nup85_seqs['nup85'], chain_id='A')
nup85_atomic0 = nup85_mol0.add_structure(nup85_pdb_fn, chain_id='A',offset=0,ca_only = True)


nup85_mol0.add_representation(nup85_mol0.get_atomic_residues(),
	resolutions=[1],
	density_residues_per_component=1,
	density_prefix=out_path+'/nup85_0_gmm',
	density_force_compute=True,
	density_voxel_size=10.0
	)

bskt_mols.append(nup85_mol0)


#-------------------------------------------
hier = s.build()
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
output_objects = []

#-------------------------------------------
for mol in bskt_mols:
	dof.create_rigid_body(mol.get_atomic_residues())

#-------------------------------------------
if include_gaussian_em_restraint:

	densities = IMP.atom.Selection(hier,
		molecule="nup85",
		copy_index = 0,
		representation_type=IMP.atom.DENSITIES).get_selected_particles()
	
	gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities, gmmfilename,
		scale_target_to_mass=True,
		slope=0.1)
	gem.set_label("EM")
	gem.add_to_model()
	gem.set_weight(weight=1.0)
	output_objects.append(gem)

	t0 = gem.evaluate()
	mdl.update()

#-------------------------------------------
gmm_com = [0, 0, 0]
BBox_Size=500
BBx1,BBy1,BBz1 = gmm_com[0]-BBox_Size/2, gmm_com[1]-BBox_Size/2, gmm_com[2]-BBox_Size/2
BBx2,BBy2,BBz2 = gmm_com[0]+BBox_Size/2, gmm_com[1]+BBox_Size/2, gmm_com[2]+BBox_Size/2
IMP.pmi.tools.shuffle_configuration(bskt_mols,
	bounding_box = ((BBx1,BBy1,BBz1),(BBx2,BBy2,BBz2)),
	max_translation=500,
	)
mdl.update()

#-------------------------------------------
dof.optimize_flexible_beads(50)
rex1 = IMP.pmi.macros.ReplicaExchange(
	mdl,
	root_hier=hier,
	monte_carlo_sample_objects=dof.get_movers(),
	global_output_directory=out_path+'/REX1/',
	output_objects=output_objects,
	write_initial_rmf=True,
	monte_carlo_steps=20,
	number_of_best_scoring_models=0,
	number_of_frames=100,
	score_moved=True,
	replica_exchange_minimum_temperature=1.0,
	replica_exchange_maximum_temperature=4.0)

rex1.execute_macro()
