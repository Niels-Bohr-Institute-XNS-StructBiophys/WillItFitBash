import os
import fnmatch
import sys
import subprocess
from shutil import copyfile, rmtree
from ntpath import basename, splitext

# compilation of C-code:
# cd C/ ; gcc WillItFit.c -O3 -fopenmp -g -Wall -o WillItFit.com -lm -lgsl -lgslcblas ; cd ..
#
# run python script:
# python WIFBash.py


InitialDirectory = os.path.abspath(os.path.dirname(sys.argv[0]))
os.chdir(InitialDirectory)

print 'InitialDirectory = ' + InitialDirectory
print ''



# Q-range in 1/\AA applied for all files
MinQStr = '0.0'
MaxQStr = '1.0'

# default arguments for FittingRoutine
#
# FittingRoutineArgument1 = 50
# FittingRoutineArgument2 = 10
# FittingRoutineArgument3 = 32
#
# 0 ComputeModel	-
#
# 1 LevenbergMarquardt	FittingRoutineArgument1
#			(MaxIterations in LM)
#
# 2 GridsearchLM	FittingRoutineArgument1		FittingRoutineArgument2
#			(MaxIterations in LM)		(NumberOfCycles in Gridsearch)
#
# 3 BFGS		FittingRoutineArgument1
#			(MaxIterations in BFGS)
#
# 4 GridsearchBFGS	FittingRoutineArgument1		FittingRoutineArgument2
#			(MaxIterations in BFGS)		(NumberOfCycles in Gridsearch)
#
# 5 Swarm		FittingRoutineArgument1		FittingRoutineArgument2		FittingRoutineArgument3
#			???				???				???
#
# 6 Genetic		-				FittingRoutineArgument2		FittingRoutineArgument3
#							???				???
#
FittingRoutine = 1
FittingRoutineArgument1 = 15
FittingRoutineArgument2 = 3	# only in use for FittingRoutine == 2, 4, 5, 6 
FittingRoutineArgument3 = 0	# only in use for FittingRoutine == 5, 6

ResolutionFile = 'N/A'
NumberOfSmearingFolds = 0

PrintCovarianceMatrix = 1 # boolean flag 0 or 1
ChiSquareFractile = 0.0

CMD = 1 # do not change
ReadAtomsAsResidues = 1 # =0 for residue-based reading of PDB and calculations, !=0 for all-atom-based (each atom is treated as a single residue)

# logging / reporting
#
# WriteLog  < 0 output printed to terminal
# WriteLog == 0 nothing logged
# WriteLog  > 0 output printed to logfile.log in results -folder
#
# abs(WriteLog) == 0 -> no output at all
#
# abs(WriteLog) == 1 -> most important output written, includes:
#			-parameter file
#			-excerpt of spectra, ProteinStructure
#			-ChiSquare and parameters at the fitting iteration steps
#
# abs(WriteLog) == 2 -> full output written, includes in addition to 1:
#			full output of Spectra, ProteinStructure
#			output from fitting algorithm (currently Levenberg-Marquardt and Compute Model supported)
#
# default is -1

WriteLog = 1


ExecutableFile = os.path.join( InitialDirectory, "C/WillItFit.com")
print('ExecutableFile = ' + ExecutableFile)


# One Data/ directory with different PDBs but same ONE *.dat, and *.par.
# Use for each job an individual *.card file, otherwise the output will be overwritten
DataDirectory = InitialDirectory + "/Data/"
print('DataDirectory = ' + DataDirectory)

SampleFile = DataDirectory + 'GHR_MSP1D1POPC.dat'
print('SampleFile = ' + SampleFile)
if not SampleFile:
	print('Sample file ' + SampleFile + ' does not exist')
	sys.exit("Exit")


ParameterFile = DataDirectory + 'best_fitAB.par'
print('ParameterFile = ' + ParameterFile)
if not ParameterFile:
	print('Parameter file ' + ParameterFile + ' does not exist')
	sys.exit("Exit")


CardTmpFile = DataDirectory + 'GHR_ND.card'
print('CardTmpFile = ' + CardTmpFile)
if not CardTmpFile:
	print('Card template file ' + CardTmpFile + ' does not exist')
	sys.exit("Exit")


# Sort PDB filelist
for PDBFile in sorted(fnmatch.filter(os.listdir( DataDirectory ), '*.pdb')):

	print ''
	print ''

	# create unique Card-file for each PDB
	CardFile = DataDirectory + splitext(PDBFile)[0] + '-' + basename(CardTmpFile)
	print CardFile
	copyfile( CardTmpFile, CardFile)

	PDBFile = DataDirectory + PDBFile
	print PDBFile

	# remove old *-results folders if exist
	print os.path.join( InitialDirectory, CardFile) + '-results'
	OldResultsFolder = os.path.join( InitialDirectory, CardFile) + '-results'
	if os.path.isdir( OldResultsFolder ) :
		rmtree( OldResultsFolder )

	# assemble WIFBash call string array
	ProcessToCall = []
	ProcessToCall.append(ExecutableFile)

	ProcessToCall.append('-c=%s' % CardFile)
	ProcessToCall.append('-s=%s' % SampleFile)
	ProcessToCall.append('-p=%s' % ParameterFile)
	ProcessToCall.append('-d=%s' % PDBFile)
	ProcessToCall.append('-a=%d' % ReadAtomsAsResidues)

	ProcessToCall.append('-n=%s' % MinQStr)
	ProcessToCall.append('-x=%s' % MaxQStr)

	ProcessToCall.append('-r=%d' % FittingRoutine)
	ProcessToCall.append('-i=%d' % FittingRoutineArgument1)
	ProcessToCall.append('-j=%d' % FittingRoutineArgument2)
	ProcessToCall.append('-k=%d' % FittingRoutineArgument3)

	ProcessToCall.append('-e=%s' % ResolutionFile)
	ProcessToCall.append('-t=%d' % NumberOfSmearingFolds)

	ProcessToCall.append('-o=%d' % PrintCovarianceMatrix)
	ProcessToCall.append('-h=%s' % ChiSquareFractile)

	ProcessToCall.append('-z=%d' % CMD)
	ProcessToCall.append('-l=%d' % WriteLog)

	print ""
	print ProcessToCall
	print ""

	# run WIFBash exec
	Process = subprocess.Popen(ProcessToCall)
	out, err = Process.communicate()
	print out, err

