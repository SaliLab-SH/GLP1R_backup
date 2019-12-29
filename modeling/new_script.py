"""This script is another script for modeling GLP1R homodimer
"""

from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.dof
import IMP.atom
import IMP.core
import IMP.algebra

import os
import sys


##### Mol Setup #####
def add_pdb_rep(mol,chain,unstructured_bead_size,clr,prot,dens):
    print(prot,chain)
    atomic = mol.add_structure(datadirectory+prot+'.pdb',chain_id=chain,offset=0)
    mol.add_representation(atomic, resolutions=[1,10],color = clr)
    print("added", prot, chain)
    if (len(mol[:]-atomic)>0):
        mol.add_representation(mol[:]-atomic,resolutions=[unstructured_bead_size],color=clr,setup_particles_as_densities=False,density_force_compute=False)
    return mol

def add_disorder_rep(mol,chain,unstructured_bead_size,clr):
    mol.add_representation(mol[:],resolutions=[unstructured_bead_size],color=clr,setup_particles_as_densities=False,density_force_compute=False)
    return mol




test_mode = False
if '--test' in sys.argv: test_mode = True

num_frames = 50000

##### Preparation #####
# Data files
directory = '/home/hqcao/Documents/PBC_project/GLP1R/Data/'
# xl_data = directory + 'xlinks.txt'
xl_data = directory + 'xl_selected_1.txt'

# Restraint weights
xl_weight = 10.0

# Bead Size
bead_size = 20

# Topology file
topology_file = directory + 'topology.txt'

##### Main #####
# Build the model
m = IMP.Model()

# Read the topology file
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro and add a state from topology file
bs = IMP.pmi.macros.BuildSystem(m)
bs.add_state(t)

# Execute the macro and get the return for root hierarchy and degrees of freedom
root_hier, dof = bs.execute_macro()
# print("Flexible beads: ", dof.get_flexible_beads())
# print("Movers: ", dof.get_movers())
# print("Rigid bodies: ", dof.get_rigid_bodies())

# sys.exit("Test DOF")


# Get the list of molecules
molecules = t.get_components()




##### Restraints #####
# Create the output list to collect the restraints
# All the restraints must be add_to_model()
output_objects = []

# Excluede Volume Restraint
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = [root_hier], resolution = 10)
ev.add_to_model()
output_objects.append(ev)

# Crosslinking Restraints
xl = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xl.set_protein1_key("prot1")
xl.set_protein2_key("prot2")
xl.set_residue1_key("res1")
xl.set_residue2_key("res2")

xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb.create_set_from_file(file_name = xl_data, converter = xl)
xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier = root_hier,
        CrossLinkDataBase = xldb,
        length = 21,
        resolution = 1.0,
        slope = 0.02,
        weight = xl_weight)
xlr.add_to_model()
output_objects.append(xlr)

##### Sampling #####
# Randomize
IMP.pmi.tools.shuffle_configuration(root_hier, max_translation = 50)
# Optimize for relaxing large connectivity
dof.optimize_flexible_beads(500)
# Run REMC sampling
rex = IMP.pmi.macros.ReplicaExchange0(m,
        root_hier = root_hier,
        crosslink_restraints = [xlr],
        monte_carlo_sample_objects = dof.get_movers(),
        global_output_directory = 'output/',
        output_objects = output_objects,
        monte_carlo_steps = 10,
        number_of_best_scoring_models = 100,
        number_of_frames = num_frames,
        test_mode = test_mode)
rex.execute_macro()

# output_xlr = open("output_xlr.txt","w")
# output_xlr.write(xlr.evaluate())
print(xlr.evaluate())
