import Bio.PDB
import sys
import string
import os
import argparse
import timeit
import logging
import re
from macrocomplex_functions import *
# alias chimera="~/.local/UCSF-Chimera64-1.13.1/bin/chimera"

start = timeit.default_timer()

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

parser.add_argument('-p', '--pdb_iterations',		
					dest = "pdb_iterations",
					action = "store_true",
					default = False,
					help = "The program create a new pdb file every time that one iteration finishes.")

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

if arguments.outdir == None:		#Check output directory. If not, creates one
	arguments.outdir = arguments.indir + "_out"	
else:
	pass
if not os.path.exists(arguments.outdir):		# Checking if the OUTPUT directory (created by us or provided by the user) already exists
	os.mkdir(arguments.outdir)		# If not, create it
	arguments.outdir = os.path.abspath(arguments.outdir)
else:
	arguments.outdir = os.path.abspath(arguments.outdir)

os.chdir(arguments.outdir)		##changes the current directory to OUTPUT directory

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

file_path = arguments.indir + "/" + files[0]		#creates path of the first file, which is going to be the reference structure
pdb_parser = Bio.PDB.PDBParser(QUIET = True)	#creation of PDBParser object
ref_structure = pdb_parser.get_structure("reference", file_path)	#creation of the reference structure with the first file, necessary to call the funtion
logging.info("The initial complex has %d chains and are the following:" % (ref_structure[0].__len__()))
for ID in [chain.get_id() for chain in ref_structure[0].get_chains()]:		#loops through all chains of ref_structure
	logging.info("Chain %s", ID)		#prints the ID

ref_structure = MacrocomplexBuilder(ref_structure = ref_structure, files_list = files, it = 0, not_added = 0, command_arguments = arguments)	#calling the iterative function

### MACROCOMPLEX BUILDING PROCESS FINISHED ###
if len(list(ref_structure[0].get_atoms())) > 99999 or len(list(ref_structure[0].get_chains())) > 62:		#checks that the structure has has less atoms than the maximum for a PDB, 99,999
	io = Bio.PDB.MMCIFIO()								#creates the MMCIFIO object, that can contain more than 99,999 atom coordinates
	io.set_structure(ref_structure[0])					#sets the reference structure object to be written in a MMCIF file
	io.save("macrocomplex.cif")							#saves the file in output directory
	logging.info("Output files %s saved in %s" %("macrocomplex.cif and macrocomplex.log",os.path.abspath(arguments.outdir)))
else: 													#checks that the structure has more than 99,999 atoms
	io = Bio.PDB.PDBIO()								#creates the PDBIO object
	io.set_structure(ref_structure[0])					#sets the reference structure object to be written in a PDB file
	io.save("macrocomplex.pdb")							#the whole macrocomplex gets saved in "macrocomplex.pdb"
	logging.info("Output files %s saved in %s" %("macrocomplex.pdb and macrocomplex.log",os.path.abspath(arguments.outdir)))
	
stop = timeit.default_timer()
logging.info("The program has finished running! It took %f seconds" % (stop - start))
