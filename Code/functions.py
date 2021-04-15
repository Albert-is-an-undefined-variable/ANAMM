import Bio.PDB
import sys
import string
import os
import argparse
import logging
import re
from modeller import environ, selection
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions


def get_key_atom(chain):
	"""This function retrieves the key atom, CA in case of proteins and C4' in case of nucleic acids, to do the superimposition and also returns a
	variable indicating the kind of molecule that that chain is: either DNA, RNA or PROTEIN
	"""
	### Declaring and creating new variables ###
	nucleic_acids = ['DA','DT','DC','DG','DI','A','U','C','G','I']		#creating a list with all possible nucleic acids letters
	RNA = ['A','U','C','G','I']											#creating a list with all possible RNA letters
	DNA = ['DA','DT','DC','DG','DI']									#creating a list with all possible DNA letters
	atoms = []
	### Loops through all residues of the chain ###
	for res in chain:
		res_name = res.get_resname()[0:3].strip()		#get the name of the residue (with no spaces)
		## Appends the CA atoms and sets the molecule type of protein ##
		if res.get_id()[0] == " " and res_name not in nucleic_acids:		#checks whether the residue is not a HETATM or nucleic acid
			if 'CA' not in res:		#checks whether the residue has CA atoms
				logging.warning("This protein residue %d %s does not have CA atom" % (res.get_id()[1], res_name))
			else:
				atoms.append(res['CA'])		#append CA atoms to the list of sample atoms
				molecule = 'PROTEIN'		#set the molecule type to protein
		## Append the C4 atoms and sets the molecule type of DNA or RNA ##
		elif res.get_id()[0] == " " and res_name in nucleic_acids:		#checks whether the residue is a nucleic acid and not HETATM
			if res_name in DNA:			#checks whether the residue is a DNA nucleotide
				molecule = 'DNA'		#set the molecule type to DNA
			elif res_name in RNA:		#checks whether the residue is a RNA nucleotide
				molecule = 'RNA'		#set the molecule type to RNA
			atoms.append(res['C4\''])	#append C4' atoms to the list of atoms
	return(atoms, molecule)		# Return all key atoms list and the type of molecule to which they belong

def chain_character_ID(ID_list, ID):
	"""
	The function outputs new characters IDs to a chain that will be included into the macrocomplex. The function requires a list with the IDs of the already existing
	chains to ensure that it does not generate an IDs that is already being used.
	"""
	characters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')	#characters for IDs

	if len(ID_list) < 62:				#if the ID list of all the chains is smaller than 62, the new chains will receive a single character ID
		if ID not in ID_list:
			return ID
		elif ID in ID_list:
			for i in characters:		#loop through characters
				if characters[i] not in ID_list:
					return characters[i]
				else:								#keep looping if the ID is already found in the ID list
					continue

	elif len(ID_list) >= 62:			#if the ID list of all the chains is bigger than 62, the new chains will receive a two-character ID
		for char1, char2 in characters:
			ID = char1 + char2	#two-character ID
			if ID not in ID_list:
				return ID
			else:				#keep looping if the ID is already found in the ID list
				continue

def optimize(pdb, pdb_path):
    print(1, pdb_path)
    # Environ data
    env = environ(0)
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    code = pdb.split('.')[0]
    mdl = complete_pdb(env, pdb)
    mdl.write(file=code+'.ini')

    # Select all atoms:
    atmsel = selection(mdl)

    # Generate the restraints:
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
    mdl.restraints.write(file=code+'.rsr')

    mpdf_prior = atmsel.energy()

    # Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='REPORT')
    md = molecular_dynamics(output='REPORT')

    # Open a file to get basic stats on each optimization
    trcfil = open(code+'.D00000001', 'w')

    # Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20, actions=actions.trace(5, trcfil))
    # Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
    # 10 steps during the run, and write stats every 10 steps
    md.optimize(atmsel, temperature=300, max_iterations=50,
                actions=[actions.write_structure(10, code+'.D9999%04d.pdb'),
                         actions.trace(10, trcfil)])
    # Finish off with some more CG, and write stats every 5 steps
    cg.optimize(atmsel, max_iterations=20,
                actions=[actions.trace(5, trcfil)])

    mpdf_after = atmsel.energy()

    mdl.write(file=os.path.join(pdb_path, 'optimized.pdb'))
    return (mpdf_prior, mpdf_after)


def superimposition(ref_structure, sample_structure, rmsd_threshold):
	"""This function, given a reference and a sample structure does the superimposition of every combination of pairs of chains and calculates the RMSD.
	It returns a dictionary with a tuple of the reference and sample chain as a tuple and the superimposer instance resulting from those two chains, as
	well as two variables, the ID of the chain with the smallest RMSD when superimposing with the reference structure and the RMSD itself
	"""
	### Saving arguments passed on to the function ###
	ref_model = ref_structure[0]			#retrieves the first and only available model of reference structure
	sample_model = sample_structure[0]		#retrieves the first and only available model of the sample structure
	### Initializing and declaring variables ###
	best_sample_chain_ID = best_ref_chain_ID = ""
	best_RMSD = 0					#variable for the lowest RMSD
	prev_RMSD = True				#variable to know we are in the first combination of pairs of chains
	superimposed_chains = False		#variable that indicates the presence of a superimposed chain (True if there is superimposed chain)
	all_superimpositions = {}		#start the dictionary that will contain all superimposition instances
	### Superimposition of every combination of pairs of chains between the reference and the sample structures ###
	## loops through all chains in the reference model ##
	for ref_chain in ref_model:
		logging.info("Processing reference chain %s", ref_chain.id)
		ref_atoms, ref_molecule = get_key_atom(ref_chain)		#Retrieves all key atoms (CA or C4') and molecule type of the sample
		## loops through all chains in the sample model ##
		for sample_chain in sample_model:
			logging.info("Processing sample chain %s", sample_chain.id)
			sample_atoms, sample_molecule = get_key_atom(sample_chain)		#Retrieves all key atoms (CA or C4') and molecule type of the sample
			if ref_molecule != sample_molecule:				#checks that the molecular types of ref chain and sample chain are the same
				logging.warning("Cannot superimpose. Reference chain %s is %s and sample chain %s is %s" %(ref_chain.get_id(), ref_molecule, sample_chain.get_id(), sample_molecule))
			elif len(ref_atoms) != len(sample_atoms):		#checks that the length of ref_atoms and sample_atoms is the same
				logging.warning("Cannot superimpose. The number of atoms of the reference chain %s is %d and the number of atoms of the sample chain %s is %d", ref_chain.get_id(), len(ref_atoms), sample_chain.get_id(), len(sample_atoms))
			## Make the superimposition between reference and sample chain ##
			else:		#everything is fine, same type of molecule, same length of atom lists
				super_imposer = Bio.PDB.Superimposer()				#creates superimposer instance
				super_imposer.set_atoms(ref_atoms, sample_atoms)	#creates ROTATION and TRANSLATION matrices from lists of atoms to align
				RMSD = super_imposer.rms 							#retrieves RMSD
				if RMSD > rmsd_threshold:
					logging.info("The RMSD between chain %s of the reference and chain %s of the sample is %f", ref_chain.id, sample_chain.id, RMSD)
					continue
				if prev_RMSD is True or RMSD < prev_RMSD:			#checks that the RMSD of this combination is smaller than the previous one
					best_sample_chain_ID = sample_chain.id
					best_ref_chain_ID = ref_chain.id 				#with this condition, the superimposer instance and other important
					best_RMSD = RMSD 								#information pertaining to the superimposition with the smallest
					prev_RMSD = RMSD 								#RMSD will be saved
				all_superimpositions[(ref_chain.id, sample_chain.id)] = super_imposer		#saving ALL superimposer instances in a dictionary
				superimposed_chains = True 							# The superimposition has been made
				logging.info("The RMSD between chain %s of the reference and chain %s of the sample is %f", ref_chain.id, sample_chain.id, RMSD)
	### checks that there has been, at least, one superimposition ###
	if superimposed_chains is True:
		all_superimpositions = sorted(all_superimpositions.items(), key=lambda k:k[1].rms)		#sorting by the lowest RMSD and saving to a list
		logging.info("The combination of chains with the lowest RMSD is ref chain %s and sample chain %s with an RMSD of %f", best_ref_chain_ID, best_sample_chain_ID, best_RMSD)
	return(all_superimpositions, superimposed_chains, best_RMSD)

def write_pdb(reference_structure):
		"""Writes a pdb file based of the specified structure"""
		io = Bio.PDB.PDBIO()
		io.set_structure(reference_structure[0])					#reference structure object written in a PDB file
		io.save("i_model.pdb")

# def stoichiometry(ref_structure, files_list, not_added, command_arguments):
# 	### Checks if the current macrocomplex satisfies the desired number of chains or just stops at iteration 150 ###
# 	chains = ref_structure[0].__len__()
# 	sto = command_arguments.stoichiometry
# 	n = not_added
# 	if chains == sto or n > len(files_list):
# 		logging.info("The whole macrocomplex has been successfully build")
# 		logging.info("The final complex has %d chains" % chains)
# 		logging.info("We have arrived to iteration %d" %(i))
# 		return 	ref_structure			#END OF THE RECURSIVE FUNCTION



def MacrocomplexBuilder(ref_structure, files_list, it, not_added, command_arguments):
	"""This recursive function superimposes the most similar chain of a binary interaction PDB file with a reference structure and adds the transformed chain to the building complex
	"""

	### Saving arguments passed on to the function ###
	i = it 															#number of iterations
	n = not_added													#number of files that have been parsed but no chain has been added
	sto = command_arguments.stoichiometry							#number of chains
	clashes_threshold = command_arguments.clashes 					#clashes threshold
	RMSD_threshold = command_arguments.rmsd_threshold 				#RMSD threshold
	indir = command_arguments.indir 								#input directory relative path
	outdir = command_arguments.outdir 								#output directory relative path

	chains = ref_structure[0].__len__()
	### Prints the current iteration and number of chains of the current complex ###
	# ref_structure = stoichiometry(ref_structure, files_list, not_added, command_arguments)

	logging.info("This is the iteration #%d of the recursive function" % i )
	logging.info("The complex has %d chains at this point" % chains)

	# ### Checks if the current macrocomplex satisfies the desired number of chains or just stops at iteration 150 ###

	if chains == sto or n > len(files_list):
		logging.info("The whole macrocomplex has been successfully build")
		logging.info("The final complex has %d chains" % chains)
		logging.info("We have arrived to iteration %d" %(i))
		return 	ref_structure			#END OF THE RECURSIVE FUNCTION

	#Select the first element of the list of files to analyze
	(sample_structure, sample_model) = process_file(files_list[0], indir)

	### Calling the superimposition function to obtain the superimposition of every combination of pairs of chains between the reference and sample structures
	all_superimpositions, superimposed_chains, best_RMSD = superimposition(ref_structure, sample_structure, RMSD_threshold)

	### There are no superimposed chains or RMSD is above the threshold --> Call again the recursive function ###
	if superimposed_chains is False or best_RMSD > RMSD_threshold:		#if condition is met, there are no superimposed chains, or the RMSD is not small enough to be considered
		file = files_list.pop(0)		#substracts the current file
		files_list.append(file)			#and adds it at the end of the list of files
		i += 1							#calling again the recursive function to analyze the next file
		n += 1
		return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)	#call again the iterative function, j does not change
	### There are superimposed chains ###
	else:
		## Loops through the superimposition dictionary, obtaining the superimposition instances and the reference and sample IDs ##
		for chains, sup in all_superimpositions:
			logging.info("We are processing the superimposition of ref chain %s with sample chain %s with an RMSD of %f" % (chains[0],chains[1], sup.rms))
			if sup.rms > RMSD_threshold:			#Checks that the superimposition has an RMSD above the threshold
				logging.info("This superimposition of ref chain %s with sample chain %s has an RMSD bigger than the threshold, therefore it is skipped" % (chains[0],chains[1]))
				continue							#if not, skip that superimposition
			sup.apply(sample_model.get_atoms())		#applies ROTATION and TRANSLATION matrices to all the atoms in the sample model
			## Gets the sample chain that was not superimposed with the reference chain --> putative chain to add ##
			chain_to_add = [chain for chain in sample_model.get_chains() if chain.get_id() != chains[0]][0]
			present_chain = False		#this variable indicates whether the chain to add is present on the building complex or not: False => not present, True => present
			sample_atoms, sample_molecule = get_key_atom(chain_to_add)		#retrieves all key atoms (CA or C4') and molecule type of chain_to_add
			logging.info("Putative chain to add is %s" % chain_to_add.id)
			## Loops through all the chains from the reference structure ##
			for chain in ref_structure[0].get_chains():
				ref_atoms, ref_molecule = get_key_atom(chain)		#retrieves all key atoms (CA or C4') and molecule type of the reference present chain
				## Makes a Neighbor Search to look for clashes between the chain to add and the chains from the reference structure ##
				Neighbor = Bio.PDB.NeighborSearch(ref_atoms)			#creates an instance of class NeighborSearch, given a list of reference atoms
				clashes = []		#declares a list that will contain all the atoms that clash between the reference and sample chains
				for atom in sample_atoms:								#loops through the list of atoms of chain_to_add
					atoms_clashed = Neighbor.search(atom.coord,5)		#produces a Neighbor search that returns all atoms/residues/chains/models/structures that have at least one atom within radius of center.
					if len(atoms_clashed) > 0:				#if there are clashes
						clashes.extend(atoms_clashed)		#adds the atoms list to the list of clashes
				if len(clashes) > clashes_threshold:		#checks that the number of total clashes is above the threshold
					present_chain = True					#then, chain_to_add is considered a chain already present in the complex
					logging.info("The number of clashes between the chain to add %s and reference chain %s is %d, therefore the chain is skipped" % (chain_to_add.id, chain.id,len(clashes)))
					break 									#skips continuing through the loop, as it already clashes with one reference chain
				## Checks that the number of total clashes is under the threshold ##
				elif len(clashes) <= clashes_threshold:
					logging.info("The number of clashes between the chain to add %s and reference chain %s is %d, it is under the threshold" % (chain_to_add.id, chain.id,len(clashes)))
					continue								#continue the loops, as we must ensure that chain_to_add does not clash with ANY reference chain
			## Rotated chain to add is not a chain already in the building macrocomplex structure, then adds it, with its original ID or with a new one ##
			if present_chain is False:
				logging.info("Chain %s superimposed with chain %s yields rotated chain %s which is not in the complex" %(chains[0],chains[1],chain_to_add.id))
				chain_ids = [chain.id for chain in ref_structure[0].get_chains()]	#list containing IDs of all chains present in reference structure
				ID = chain_character_ID(chain_ids, chain_to_add.id)
				chain_to_add.id = ID
				ref_structure[0].add(chain_to_add)	#adds chain_to_add to the building macrocomplex structure
				logging.info("Added Chain %s" % ID)
				file = files_list.pop(0)		#substracts the first file of the files list
				files_list.append(file)			#adds the file at the end of the files list
				i += 1							#adds one to the iteration variable
				n = 0
				#this is what makes the function recursive, it calls itself on the return, executing the whole function again and again until certain condition is met
				return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)
	### Once the current file has been analyzed it is substracted and appended at the end of the files list ###
	file = files_list.pop(0)		#substracts the first file of the files list
	files_list.append(file)			#adds the file at the end of the files list
	i += 1							#adds one to the iteration variable
	n += 1
	#this is what makes the function recursive, it calls itself on the return, executing the whole function again and again until certain condition is met
	return MacrocomplexBuilder(ref_structure = ref_structure, files_list = files_list, it = i, not_added = n, command_arguments = command_arguments)

def process_file (file, input_dir):
	"""Function that requires a file and its path in order to process it. The processing consists in creating a PDBParser object and obtaining the structure and model from there.
	"""
	logging.info("We are processing the file %s" % (file))
	file_path = input_dir + "/" + file 			#path of the file
	pdb_parser = Bio.PDB.PDBParser(QUIET = True)		#obtain PDBParser object
	file_structure = pdb_parser.get_structure("file", file_path)		#structure
	file_model = file_structure[0]			#model
	return (file_structure, file_model)
