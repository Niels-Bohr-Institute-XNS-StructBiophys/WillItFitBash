// basically only Parameter, Data, ProteinStructure, UserDefinedStructure are provided 
void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters, struct Dataset * Data, int NumberOfSpectra, char *cardfilename, struct Protein ProteinStructure, struct UserDefined UserDefinedStructure, char *SampleFilename, char *ResultsDirectory)
{
	/// Declarations

	// Variables describing the output file
	FILE *fp ;
	char Filename[256] ;
	sprintf(Filename, "%sResults.wif", ResultsDirectory) ;

	/// I/O
	// Create the output file
	fp = fopen(Filename, "w+") ;

	// Print filenames to file
	fprintf( fp, ".card-file:\n") ;
	fprintf( fp, "\t%s \n", cardfilename) ;

	fprintf( fp, "\n") ;
	fprintf( fp, "Datafiles:\n") ;

	for( int i = 0; i < NumberOfSpectra; ++i) { fprintf( fp, "\t%s \n", Data[i].Filename) ; }
	fprintf( fp, "\n") ;

	// Print location of sample file
	fprintf( fp, "Sample info-file:\n") ;
	fprintf( fp, "\t%s \n", SampleFilename) ;
	fprintf( fp, "\n") ;

	// Print fit quality to file
	fprintf( fp, "Final chisq = %g \n", ChiSquare) ;
	fprintf( fp, "\n") ;

	// Print range of q to file
	fprintf( fp, "Lower limit on q = %g \n", QMin) ;
	fprintf( fp, "Upper limit on q = %g \n", QMax) ;
	fprintf( fp, "\n") ;

	// Print parameters and properties of parameters to file
	fprintf( fp, "Parameters:\n\n") ;
	fprintf( fp, "\t     Value             Error             Name\n") ;
	for( int i = 0; i < NumberOfParameters; ++i)
	{
		if ( Parameters[i].iParameter == true ) { fprintf( fp, "\t%2d   %-15g   %-15g   %s\n", i, Parameters[i].Value, Parameters[i].Error, Parameters[i].Name) ; }
		else { fprintf( fp, "\t%2d   %-15g   Fixed             %s\n", i, Parameters[i].Value, Parameters[i].Name) ; }
	}
	fprintf( fp, "\n\n") ;

	// Print info on protein file
	fprintf( fp, "Location of inital .pdb-file: %s\n", ProteinStructure.PDBFileLocation);
	fprintf( fp, "\n");

	fprintf( fp, "Atom and residue counts:\n") ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Atoms                    : %d\n", ProteinStructure.NumberOfAtoms) ;
	fprintf( fp, "\tNumber of Atoms (Modification %s) : %d\n", ProteinStructure.ModificationName, ProteinStructure.NumberOfModificationAtoms) ;
	fprintf( fp, "\tNumber of Residues                 : %d\n", ProteinStructure.NumberOfResidues) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of H  : %d\n", ProteinStructure.NumberOfHAtoms) ;
	fprintf( fp, "\tNumber of D  : %d\n", ProteinStructure.NumberOfDAtoms) ;
	fprintf( fp, "\tNumber of C  : %d\n", ProteinStructure.NumberOfCAtoms) ;
	fprintf( fp, "\tNumber of N  : %d\n", ProteinStructure.NumberOfNAtoms) ;
	fprintf( fp, "\tNumber of O  : %d\n", ProteinStructure.NumberOfOAtoms) ;
	fprintf( fp, "\tNumber of P  : %d\n", ProteinStructure.NumberOfPAtoms) ;
	fprintf( fp, "\tNumber of S  : %d\n", ProteinStructure.NumberOfSAtoms) ;
	// fprintf( fp, "Number of I: %d\n", ProteinStructure.NumberOfIAtoms) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Zn : %d\n", ProteinStructure.NumberOfZNAtoms) ;
	fprintf( fp, "\tNumber of Cl : %d\n", ProteinStructure.NumberOfCLAtoms) ;
	fprintf( fp, "\tNumber of Na : %d\n", ProteinStructure.NumberOfNAAtoms) ;
	fprintf( fp, "\tNumber of Ca : %d\n", ProteinStructure.NumberOfCAAtoms) ;
	fprintf( fp, "\tNumber of Fe : %d\n", ProteinStructure.NumberOfFEAtoms) ;
	fprintf( fp, "\n") ;
	fprintf( fp, "\tNumber of Q  : %d\n", ProteinStructure.NumberOfQAtoms) ;
	fprintf( fp, "\n\n") ;

	// Print constraints to file
	for ( int i = 0; i < NumberOfSpectra; ++i)
	{
		fprintf( fp, "Constraints for spectrum %d: \n", i);
		fprintf( fp, "\tScattering length density of headgroups = %g\n", Data[i].Constraints[0]);
		fprintf( fp, "\tScattering length density of core       = %g\n", Data[i].Constraints[1]);
		fprintf( fp, "\tScattering length density of methyl     = %g\n", Data[i].Constraints[2]);
		fprintf( fp, "\tScattering length density of belt       = %g\n", Data[i].Constraints[3]);
		fprintf( fp, "\tArea per headgroup                      = %g\n", Data[i].Constraints[9]);
		fprintf( fp, "\n");
	}

	// Print global parameters to file
	fprintf( fp, "Global parameteres:\n");
	fprintf( fp, "\tHeight of bilayer             = %g \n", Data[0].Constraints[6]);
	fprintf( fp, "\tHeight of hydrophobic bilayer = %g \n", Data[0].Constraints[7]);
	fprintf( fp, "\tHeight of methyl layer        = %g \n", Data[0].Constraints[8]);
	fprintf( fp, "\tConcentration                 = %g \n", Data[0].Constraints[12]);
	fprintf( fp, "\tMol. volume of dry lipid head = %g \n", Data[0].Constraints[13]);
	fprintf( fp, "\tMol. volume of lipid tail     = %g \n", Data[0].Constraints[14]);
	fprintf( fp, "\tMol. volume of methyl         = %g \n", Data[0].Constraints[15]);
	fprintf( fp, "\tMol. volume of belt           = %g \n", Data[0].Constraints[16]);
	fprintf( fp, "\n");
	fprintf( fp, "\tNumber density       = %g \n", Data[0].Constraints[21]);
	fprintf( fp, "\tMajor semiaxis       = %g \n", Data[0].Constraints[23]);
	fprintf( fp, "\tMinor semiaxis       = %g \n", Data[0].Constraints[24]);
	fprintf( fp, "\tVolume of core       = %g \n", Data[0].Constraints[25]);
	fprintf( fp, "\tVolume of headgroups = %g \n", Data[0].Constraints[26]);
	fprintf( fp, "\tWidth of belt        = %g \n", Data[0].Constraints[27]);
	fprintf( fp, "\tVolume of methyl     = %g \n", Data[0].Constraints[28]);

	// Close file
	fclose(fp) ;


	struct Residue CurrentResidue;
	double MPTranslation[3];
	double MPRotation[3][3];
	SetTranslationVector(MPTranslation, Parameters[24].Value, Parameters[25].Value, Parameters[26].Value);
	SetRotationMatrix(MPRotation, Parameters[27].Value, Parameters[28].Value, Parameters[29].Value);

	sprintf( Filename, "%sout.pdb", ResultsDirectory) ;
	fp = fopen( Filename, "w+") ;

	for( int j = 0; j < ProteinStructure.NumberOfResidues; ++j)
	{
		CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
		Orient(&CurrentResidue.xVolume, &CurrentResidue.yVolume, &CurrentResidue.zVolume, MPRotation, MPTranslation);

		fprintf( fp, "ATOM   %4d  CA  %s   %4d    %7.3lf %7.3lf %7.3lf  1.00 18.20           N\n", j + 1, CurrentResidue.Name, j + 1, CurrentResidue.xVolume, CurrentResidue.yVolume, CurrentResidue.zVolume);
	}

	fclose(fp) ;
}
