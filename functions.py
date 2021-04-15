# Notice that we import each of the functions that we use
# this will help us to make easier to keep track of the errors
import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Superimposer
from Bio.PDB import NeighborSearch
import sys
import string
import os
import argparse
from modeller import environ, selection
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions


def fetch_backbone(chain, verbose):
    """This function retrieves the key atoms (CA for proteins, C4' for nucleic acids), to do the superimposition and also returns a
    variable indicating the kind of molecule that that chain is: either DNA, RNA or PROTEIN
    """
    # Declaring and creating new variables
    nucleic_acids = ['DA', 'DT', 'DC', 'DG', 'DI', 'A', 'U', 'C', 'G', 'I']
    RNA = ['A', 'U', 'C', 'G', 'I']
    DNA = ['DA', 'DT', 'DC', 'DG', 'DI']
    atoms = []
    # Loop through residues of the chain
    for res in chain:
        res_name = res.get_resname()[0:3].strip()  # get the name of the residue
        # Append CA atoms and set the molecule type to protein
        if res.get_id()[0] == " " and res_name not in nucleic_acids:
            if 'CA' not in res:
                if verbose:
                    print("WARNING: This protein residue %d %s does not have CA atom" %
                          (res.get_id()[1], res_name))
            else:
                atoms.append(res['CA'])
                molecule = 'PROTEIN'

        ## Append the C4 atoms and sets the molecule type of DNA or RNA ##
        elif res.get_id()[0] == " " and res_name in nucleic_acids:
            if res_name in DNA:
                molecule = 'DNA'
            elif res_name in RNA:
                molecule = 'RNA'
            atoms.append(res['C4\''])
    return(atoms, molecule)		# Return all key atoms list and the type of molecule to which they belong


def chain_character_ID(IDs, ID):
    """
    The function outputs new characters IDs to a chain that will be included into the macrocomplex. The function requires a list with the IDs of the already existing
    chains to ensure that it does not generate an IDs that is already being used.
    """
    UP = list(string.ascii_uppercase)
    LOW = list(string.ascii_lowercase)
    DIG = list(string.digits)
    # creates an alphabet containing all the possible characters that can be used as chain IDs
    alphabet = UP + LOW + DIG

    if len(IDs) < 62:
        if ID not in IDs:
            return ID
        elif ID in IDs:  # checks if the ID by default is indeed on the list
            for i in range(0, len(alphabet)):
                if alphabet[i] not in IDs:
                    return alphabet[i]
                else:  # if it is already an ID, keeps looping through the alphabet
                    continue
    elif len(IDs) >= 62:  # checks if the length of the lsit containing the IDs of all chains in the complex is greater than 62
        for char1 in alphabet:
            for char2 in alphabet:
                ID = char1 + char2
                if ID not in IDs:
                    return ID
                else:
                    continue


def optimize(pdb, pdb_path):
    """The function try to optimize the macrocomplex final model using energy profiles. It applies some general restraints using modeller package.
    """
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
    atmselection = selection(mdl)

# Generate the restraints:
    mdl.restraints.make(atmselection, restraint_type='stereo', spline_on_site=False)
    mdl.restraints.write(file=code+'.rsr')

    mpdf_prior = atmselection.energy()

# Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='REPORT')
    md = molecular_dynamics(output='REPORT')

# Open a file to get basic stats on each optimization
    trcfil = open(code+'.D00000001', 'w')

# Run CG on the all-atom selection; write stats every 5 steps
    cg.optimize(atmselection, max_iterations=20, actions=actions.trace(5, trcfil))
# Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
# 10 steps during the run, and write stats every 10 steps
    md.optimize(atmselection, temperature=300, max_iterations=50,
                actions=[actions.write_structure(10, code+'.D9999%04d.pdb'),
                         actions.trace(10, trcfil)])
# Finish off with some more CG, and write stats every 5 steps
    cg.optimize(atmselection, max_iterations=20,
                actions=[actions.trace(5, trcfil)])

    mpdf_after = atmselection.energy()

    mdl.write(file=os.path.join(pdb_path, 'optimized.pdb'))
    return (mpdf_prior, mpdf_after)


def superimposition(reference_bb, target_bb, rmsd_threshold, verbose):
    """This function, takes as arguments a sample structure and a reference structure and performs the superimposition of every combination of pairs of chains and calculates the RMSD.
    The output is a dictionary with a tuple of the reference and sample chain including the resulting superimposer instance, as
    well as two additional variables, the ID of the chain with the smallest RMSD and the RMSD value
    """
    # Saving arguments passed on to the function
    ref_model = reference_bb[0]
    target_model = target_bb[0]
    # Initializing and declaring variables
    best_sample_chain_ID = best_ref_chain_ID = ""
    best_RMSD = 0
    prev_RMSD = True
    superimposed_chains = False
    all_superimpositions = {}
    # Superimposition of every combination of pairs of chains between the reference and the sample structures
    # loop through reference model chains
    for ref_chain in ref_model:
        if verbose:
            print("Processing reference chain %s", ref_chain.id)
        ref_atoms, ref_molecule = fetch_backbone(ref_chain, verbose)
        # loop through sample model chains
        for sample_chain in target_model:
            if verbose:
                print("Processing sample chain %s", sample_chain.id)
            sample_atoms, sample_molecule = fetch_backbone(sample_chain, verbose)
            if ref_molecule != sample_molecule:  # see that the molecular types of ref chain and sample chain are the same
                if verbose:
                    print("Cannot superimpose. Reference chain %s is %s and sample chain %s is %s" % (
                        ref_chain.get_id(), ref_molecule, sample_chain.get_id(), sample_molecule))
            # ensure that the length in atoms of ref_chain and sample_chain is the same
            elif len(ref_atoms) != len(sample_atoms):
                if verbose:
                    print("Cannot superimpose. The number of atoms of the reference chain %s is %d and the number of atoms of the sample chain %s is %d",
                          ref_chain.get_id(), len(ref_atoms), sample_chain.get_id(), len(sample_atoms))
            # Make the superimposition
            else:
                super_imposer = Bio.PDB.Superimposer()
                super_imposer.set_atoms(ref_atoms, sample_atoms)
                RMSD = super_imposer.rms
                if RMSD > rmsd_threshold:
                    if verbose:
                        print("The RMSD between chain %s of the reference and chain %s of the sample is %f",
                              ref_chain.id, sample_chain.id, RMSD)
                    continue
                if prev_RMSD is True or RMSD < prev_RMSD:  # checks that the RMSD of this combination is smaller than the previous one
                    best_sample_chain_ID = sample_chain.id
                    best_ref_chain_ID = ref_chain.id
                    best_RMSD = RMSD
                    prev_RMSD = RMSD
                # saving all superimposer instances
                all_superimpositions[(ref_chain.id, sample_chain.id)] = super_imposer
                superimposed_chains = True
                if verbose:
                    print("The RMSD between chain %s of the reference and chain %s of the sample is %f",
                          ref_chain.id, sample_chain.id, RMSD)
    # ensure one superimposition
    if superimposed_chains is True:
        all_superimpositions = sorted(all_superimpositions.items(
        ), key=lambda k: k[1].rms)  # sorting by the lowest RMSD
        if verbose:
            print("The combination of chains with the lowest RMSD is ref chain %s and sample chain %s with an RMSD of %f",
                  best_ref_chain_ID, best_sample_chain_ID, best_RMSD)
    return(all_superimpositions, superimposed_chains, best_RMSD)


def write_pdb(reference_structure):
    """Writes a pdb file based of the specified structure"""
    io = Bio.PDB.PDBIO()
    # If the file has more than 99999 atoms PDBIO rises an error. Nevertheless this
    # is not a problem since the pdb file is generated.
    # We could have solved this by saiving it on a .cif file but, we have choose to
    # not do such thing in order to make easier further optimization process
    io.set_structure(reference_structure[0])  # reference structure object written in a PDB file
    io.save("created_model.pdb")


def process_file(file, input_dir, verbose):
    """Function that requires a file and its path in order to process it. The processing consists in creating a PDBParser object and obtaining the structure and model from there.
    """
    if verbose:
        print("We are processing the file %s" % (file))
    file_path = input_dir + "/" + file
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)  # PDBParser object
    file_structure = pdb_parser.get_structure("file", file_path)  # structure
    file_model = file_structure[0]  # model
    return (file_structure, file_model)


def MacrocomplexBuilder(reference_bb, files_list, it, not_added, command_arguments):
    """This recursive function superimposes the most similar chain of a binary interaction PDB file with a reference structure and adds the transformed chain to the building complex
    """

    ### Saving arguments passed on to the function ###
    i = it
    n = not_added
    st = command_arguments.stoichiometry
    clashes_threshold = command_arguments.clashes
    threshold_RMSD = command_arguments.rmsd_threshold
    indir = command_arguments.indir
    outdir = command_arguments.outdir
    verbose = command_arguments.verbose

    chains = reference_bb[0].__len__()

    # Checks if the current macrocomplex satisfies the desired number of chains or just stops

    if chains == st or n > len(files_list):
        if verbose:
            print("The whole macrocomplex has been successfully build")
            print("The final complex has %d chains" % chains)
            print("We have arrived to iteration %d" % (i))
        return reference_bb  # END

    # Select the first element of the list of files to analyze
    (target_bb, target_model) = process_file(files_list[0], indir, verbose)

    # Calling the superimposition function
    all_superimpositions, superimposed_chains, best_RMSD = superimposition(
        reference_bb, target_bb, threshold_RMSD, verbose)

    # There are no superimposed chains or RMSD is above the threshold then call the function again
    if superimposed_chains is False or best_RMSD > threshold_RMSD:
        file = files_list.pop(0)
        files_list.append(file)
        i += 1  # calling again to analyze the next file
        n += 1
        return MacrocomplexBuilder(reference_bb=reference_bb, files_list=files_list, it=i, not_added=n, command_arguments=command_arguments)
    # There are superimposed chains
    else:
        # Loop through the superimposition dictionary
        for chains, sup in all_superimpositions:
            if verbose:
                print("We are processing the superimposition of ref chain %s with sample chain %s with an RMSD of %f" % (
                    chains[0], chains[1], sup.rms))
            if sup.rms > threshold_RMSD:
                if verbose:
                    print("This superimposition of ref chain %s with sample chain %s has an RMSD bigger than the threshold, therefore it is skipped" % (
                        chains[0], chains[1]))
                continue
            # applies rotation and translation matrices to all the atoms in the sample model
            sup.apply(target_model.get_atoms())
            # Get the sample chain that was not superimposed with the reference chain and add it as a putative chain
            chain_to_add = [chain for chain in target_model.get_chains() if chain.get_id()
                            != chains[0]][0]
            present_chain = False
            sample_atoms, sample_molecule = fetch_backbone(chain_to_add, verbose)
            # Loops through all the chains from the reference structure
            for chain in reference_bb[0].get_chains():
                ref_atoms, ref_molecule = fetch_backbone(chain, verbose)
                # Makes a Neighbor Search to look for clashes between the chain to add and the chains from the reference structure
                Neighbor = Bio.PDB.NeighborSearch(ref_atoms)
                clashes = []
                for atom in sample_atoms:
                    atoms_clashed = Neighbor.search(atom.coord, 5)
                    if len(atoms_clashed) > 0:
                        clashes.extend(atoms_clashed)
                if len(clashes) > clashes_threshold:
                    present_chain = True
                    if verbose:
                        print("The number of clashes between the chain to add %s and reference chain %s is %d, therefore the chain is skipped" % (
                            chain_to_add.id, chain.id, len(clashes)))
                    break  # skips in the case that it already clashes with one reference chain
                # Checks that the number of total clashes is under the threshold
                elif len(clashes) <= clashes_threshold:
                    if verbose:
                        print("The number of clashes between the chain to add %s and reference chain %s is %d, it is under the threshold" % (
                            chain_to_add.id, chain.id, len(clashes)))
                    continue
            # the rotated chain is not a chain already in the building macrocomplex structure, then adds it
            if present_chain is False:
                if verbose:
                    print("Chain %s superimposed with chain %s yields rotated chain %s which is not in the complex" % (
                        chains[0], chains[1], chain_to_add.id))
                chain_ids = [chain.id for chain in reference_bb[0].get_chains()]
                ID = chain_character_ID(chain_ids, chain_to_add.id)
                chain_to_add.id = ID
                reference_bb[0].add(chain_to_add)
                if verbose:
                    print("Added Chain %s" % ID)
                file = files_list.pop(0)
                files_list.append(file)
                i += 1
                n = 0
                # call function again until the conditions are fulfilled
                return MacrocomplexBuilder(reference_bb=reference_bb, files_list=files_list, it=i, not_added=n, command_arguments=command_arguments)
    # the anylzed file is substracted and appended to the end of the file list
    file = files_list.pop(0)
    files_list.append(file)
    i += 1
    n += 1
    # call function again until the conditions are fulfilled
    return MacrocomplexBuilder(reference_bb=reference_bb, files_list=files_list, it=i, not_added=n, command_arguments=command_arguments)
