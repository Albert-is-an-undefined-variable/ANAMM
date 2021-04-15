import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Superimposer
from Bio.PDB import NeighborSearch
import sys
import os
import argparse
import logging
import re
from functions import *

parser = argparse.ArgumentParser(description =
	"The function of this program is to reconstruct biological macrocomplexes. It can build them using protein-DNA/RNA interactions as well as protein-protein interactions. The input is a set of binary interactions and the desired number of chains of the target complex. Moreover, you can add a RMSD_treshold and the number of pdb_iteration as extra arguments. As input the program needs the folder where the input files are, and you have to select an output directory to save the output files. If you do not select the output folder, the program automatically will create one.")


Mandatory = parser.add_argument_group('Mandatory arguments')

Mandatory.add_argument('-in', '--indir',
					dest = "indir",
					action = "store",
					required=True,
					help = "MANDATORY ARGUMENT: Directory containing all PDB files with the protein binary interactions.")

parser.add_argument('-out', '--outdir',
					dest = "outdir",
					action = "store",
					default = None,
					help = " Directory where final results will be save. The program generates a macrocomplex.pdb (file with the final structure) and macrocomplex.log (log file). By default, if the user does not indicate any output directory, the program will create one. The name of the output created by default by the program is the input folder * \"_out\".")

parser.add_argument('-v', '--verbose',
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "If set, the progression log printed in standard output file.")


parser.add_argument('-sto', '--stoichiometry',
					dest = "stoichiometry",
					action = "store",
					type = int,
					default = 100,
					help = "Argument that allows you to indicate the stoichiometry desired for the target complex, indicating the final number of chains that you macrocomplex have.")

parser.add_argument('-rt', '--rmsd_threshold',
					dest = "rmsd_threshold",
					action = "store",
					default = 0.3,
					type = float,
					help = "WARNING: Argument misuse could affect the final output. If set, the RMSD threshold for considering a superimposition as correct will take this value. If not, it will be 0.3 by default.")

parser.add_argument('-cl', '--clashes_threshold',
					dest = "clashes",
					action = "store",
					default = 30,
					type = int,
					help = "WARNING: Argument misuse could affect the final output. If set, the threshold of the number of clashes will take this value. If not, it will be 30 by default.")

parser.add_argument('-opt','--optimize',
                    action='store_true',
                    help='Whether the model will be optimized with MODELLER after building or not. Default False',
                	dest='optimize',
                    default=False)

arguments = parser.parse_args()  #saving command-line arguments

if __name__ == "__main__":

	if not arguments.indir:
		raise NameError("ERROR! YOU NEED TO INDICATE INPUT DIRECTORY. Please, use the help flag, --help, h if you need instructions!")
	else:
		if (os.path.isdir(arguments.indir)):
			files = sorted(list(filter(lambda x: x.endswith(".pdb"), os.listdir(arguments.indir))))	#keep in a list of files only the files ending with .pdb
			arguments.indir = os.path.abspath(arguments.indir)
		else:
			raise NameError("ERROR! INPUT DIRECTORY NOT FOUND. CHECK THE NAME OR THE PATH!")

	os.chdir(arguments.indir + "/../")

	if arguments.outdir == None:		#look for output directory. If not, create it
		arguments.outdir = arguments.indir + "_out"
	else:
		pass
	if not os.path.exists(arguments.outdir):		# looking if the output directory already exists
		os.mkdir(arguments.outdir)		# if not, create it
		arguments.outdir = os.path.abspath(arguments.outdir)
	else:
		arguments.outdir = os.path.abspath(arguments.outdir)

	os.chdir(arguments.outdir)		#changes the current directory to output directory

	# Initializing the LOG system

	logging.basicConfig(format = '%(levelname)s:%(message)s', filename = arguments.outdir + '/macrocomplex.log', level = logging.DEBUG)
	logging.debug('...STARTING...')		# The LOG file is "macrocomplex.log" by default

	if arguments.verbose:		# Checking if VERBOSE argument is set
		logging.getLogger().addHandler(logging.StreamHandler())		# If it is set, the contents of the LOG file is also printed in STDOUT

	if re.search("\d",files[0]):
		files = sorted(files, key=lambda x: str("".join([i for i in x if i.isdigit()])))
	else:
		files = sorted(files)

	logging.info("Parameters used are:\n  - Stoichiometry: %d\n  - RMSD threshold: %.4f\n  - Clashes threshold %d" % (arguments.stoichiometry, arguments.rmsd_threshold, arguments.clashes))

	#Use the structure of the first file as the reference
	(ref_structure, ref_model_useless) = process_file(files[0], arguments.indir) #Function returns model but we are only interested in the structure


	logging.info("The initial complex has %d chains and are the following:" % (ref_structure[0].__len__()))
	for ID in [chain.get_id() for chain in ref_structure[0].get_chains()]:		#loops through all chains of ref_structure
		logging.info("Chain %s", ID)		#prints the ID

	ref_structure = MacrocomplexBuilder(ref_structure = ref_structure, files_list = files, it = 0, not_added = 0, command_arguments = arguments)	#calling the iterative function


	write_pdb(ref_structure)

	logging.info("Output files %s saved in %s" %
	             ("macrocomplex.pdb and macrocomplex.log", os.path.abspath(arguments.outdir)))
	print(f"Finished! you can see the results in {arguments.outdir}")


	# Try to optimize the final complex
	if(arguments.optimize):
		from builder import optimize
		print('Optimizing...')
		sys.stdout = open(os.devnull, 'w')
		energies = optimize.optimize(os.path.join(arguments.outdir, "i_model.pdb")	, arguments.outdir)
		sys.stdout = sys.__stdout__
		print('Energy before optimizing: %s' % str(energies[0][0]))
		print('Energy after optimizing: %s' % str(energies[1][0]))
		print('Model completed')
