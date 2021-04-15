import Bio.PDB
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


parser.add_argument('-nc', '--num_chains',
					dest = "num_chains",
					action = "store",
					type = int,
					default = 100,
					help = "Argument that allows you to indicate the number of chains desired for the target complex.")

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

arguments = parser.parse_args()  #Saving and checking the command-line arguments


### CHECKING CORRECT USE OF ARGUMENTS  ###

if not arguments.indir:
	raise NameError("ERROR! YOU NEED TO INDICATE INPUT DIRECTORY. Please, use the help flag, --help, h if you need instructions!")
else:
	if (os.path.isdir(arguments.indir)):
		files = sorted(list(filter(lambda x: x.endswith(".pdb"), os.listdir(arguments.indir))))	#Keep in a list of files only the files ending with .pdb
		arguments.indir = os.path.abspath(arguments.indir)
	else:
		raise NameError("ERROR! INPUT DIRECTORY NOT FOUND. CHECK THE NAME OR THE PATH!")

os.chdir(arguments.indir + "/../")

if arguments.outdir == None:		#Check output directory. 
	arguments.outdir = arguments.indir + "_out"  #If not, creates one
else:
	pass
if not os.path.exists(arguments.outdir):		# Checking if the OUTPUT directory already exists
	os.mkdir(arguments.outdir)		# If not, create it
	arguments.outdir = os.path.abspath(arguments.outdir)
else:
	arguments.outdir = os.path.abspath(arguments.outdir)

os.chdir(arguments.outdir)		#makes the current directory into the output directory

### Initializing the LOG system ###

logging.basicConfig(format = '%(levelname)s:%(message)s', filename = arguments.outdir + '/macrocomplex.log', level = logging.DEBUG)
logging.debug('...STARTING...')		# The LOG file is "macrocomplex.log" by default

if arguments.verbose:		# Checking if VERBOSE argument is set
	logging.getLogger().addHandler(logging.StreamHandler())		# If it is set, the LOG file is also printed in STDOUT

if re.search("\d",files[0]):		# Checking if the files
	files = sorted(files, key=lambda x: str("".join([i for i in x if i.isdigit()])))
else:
	files = sorted(files)

logging.info("Parameters used are:\n  - Number of chains: %d\n  - RMSD threshold: %.4f\n  - Clashes threshold %d" % (arguments.num_chains, arguments.rmsd_threshold, arguments.clashes))

#first file reference

file_path = arguments.indir + "/" + files[0]		#path of the first file, the reference
pdb_parser = Bio.PDB.PDBParser(QUIET = True)	
ref_structure = pdb_parser.get_structure("reference", file_path)	#reference structure from the first file
logging.info("The initial complex has %d chains and are the following:" % (ref_structure[0].__len__()))
for ID in [chain.get_id() for chain in ref_structure[0].get_chains()]:		#loops through all chains of the reference structure
	logging.info("Chain %s", ID)		#print ID

ref_structure = MacrocomplexBuilder(ref_structure = ref_structure, files_list = files, it = 0, not_added = 0, command_arguments = arguments)	#FIRST call of the iterative function


io = Bio.PDB.PDBIO()								
io.set_structure(ref_structure[0])					#reference structure object written in a PDB file
io.save("macrocomplex.pdb")	
#save in "macrocomplex.pdb"
#logging.info("Output files %s saved in %s" %("macrocomplex.pdb and macrocomplex.log",os.path.abspath(arguments.outdir)))

#logging.info("The program has finished running! It took %f seconds" % (stop - start))


# Try to optimize the build   #

if(arguments.optimize):
	from builder import optimize
	print('Optimizing...')
	sys.stdout = open(os.devnull, 'w')
	energies = optimize.optimize(os.path.join(arguments.outdir, 'macrocomplex.pdb'), arguments.outdir)
	sys.stdout = sys.__stdout__
	print('Energy before optimizing: %s' % str(energies[0][0]))
	print('Energy after optimizing: %s' % str(energies[1][0]))
	print('Model completed')
