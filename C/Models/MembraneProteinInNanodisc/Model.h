double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{

    // Get scattering length densities
    double ScatteringLengthDensityOfCaps   = Constraints[0];
    double ScatteringLengthDensityOfCore   = Constraints[1];
    double ScatteringLengthDensityOfMethyl = Constraints[2];
    double ScatteringLengthDensityOfBelt   = Constraints[3];
    double ScatteringLengthDensityOfWater  = Constraints[4];
    
    // Get the parameters describing the geometry of the disc
    double HeightOfBelt        = Parameters[2];
    double HeightOfLipids      = Constraints[6];
    double HeightOfCore        = Constraints[7];
    double HeightOfMethyl      = Constraints[8];
    double MajorSemiAxisOfCore = Constraints[23];
    double MinorSemiAxisOfCore = Constraints[24];
    double ThicknessOfBelt     = Constraints[27];

    // Get the parameters describng the membrane protein
    double CorrectionToMPVolume = Parameters[11];
    double MPTranslation[3];
    double MPRotation[3][3];
    SetTranslationVector(MPTranslation, Parameters[24], Parameters[25], Parameters[26]);
    SetRotationMatrix(MPRotation, Parameters[27], Parameters[28], Parameters[29]);

    // Parameters describing the solution
    double Roughness;
    double Scaling;
    double Background;
    double Intensity;
    double ConcentrationOfSample = Constraints[21];
   
    int j,l,m;
    struct Residue CurrentResidue;

	// Structures for holding spherical harmonics
    double complex ** Beta  = ComplexArray(Nh + 1, Nh + 1);
    double complex ** Alpha = ComplexArray(Nh + 1, Nh + 1);

    // Compute the roughness from the parameters
    if (Contrast >= 0.0 && Contrast <= 100.0) {
        Roughness  = exp(-(q * Parameters[18] * q * Parameters[18]));
        Scaling    = Parameters[20] * (Contrast / 100.0 * Parameters[13] + (100.0 - Contrast) / 100.0 * Parameters[15]);
        Background = Contrast / 100.0 * Parameters[12] + (100.0 - Contrast) / 100.0 * Parameters[14];
    } else {
        Roughness  = exp(-(q * Parameters[8] * q * Parameters[8]));
        Scaling    = Parameters[17] * Parameters[20];
        Background = Parameters[16];
    }

    // Nanodisc
    NanodiscPDBModel(Alpha,q, ScatteringLengthDensityOfCaps, ScatteringLengthDensityOfCore, ScatteringLengthDensityOfMethyl,ScatteringLengthDensityOfBelt, ScatteringLengthDensityOfWater, MajorSemiAxisOfCore, MinorSemiAxisOfCore, ThicknessOfBelt, HeightOfLipids, HeightOfCore, HeightOfMethyl, HeightOfBelt);

    // Protein structure
    for (j = 0; j < ProteinStructure.NumberOfResidues; ++j) {
        CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
        Orient(&CurrentResidue.xVolume, &CurrentResidue.yVolume, &CurrentResidue.zVolume, MPRotation, MPTranslation);

        if (Contrast < 0.0 || Contrast > 100.0) {
            Orient(&CurrentResidue.xXRayScattering, &CurrentResidue.yXRayScattering, &CurrentResidue.zXRayScattering, MPRotation, MPTranslation);
        } else {
            Orient(&CurrentResidue.xNeutronScattering, &CurrentResidue.yNeutronScattering, &CurrentResidue.zNeutronScattering, MPRotation, MPTranslation);
        }

       AddScatteringFromResidue(Beta, q, CurrentResidue, Contrast, ScatteringLengthDensityOfCaps, ScatteringLengthDensityOfCore, ScatteringLengthDensityOfMethyl, ScatteringLengthDensityOfBelt, ScatteringLengthDensityOfWater, HeightOfLipids, HeightOfCore, HeightOfMethyl, HeightOfBelt, CorrectionToMPVolume, MajorSemiAxisOfCore, MinorSemiAxisOfCore);
    }

    // Calculate intensity
    Intensity = 0.0;

    for (l = 0; l < Nh + 1; ++l) {

        for (m = 0; m < l + 1; ++m) {
            Intensity += ((m > 0) + 1) * pow(cabs(sqrt(Roughness) * Alpha[l][m] + Beta[l][m]), 2);
        }
    }

    // Scale the sum
    Intensity *= ConcentrationOfSample;
    Intensity = Intensity * Scaling + Background;

    // Free the arrays
    FreeComplexArray(Alpha, Nh + 1, Nh + 1);
    FreeComplexArray( Beta, Nh + 1, Nh + 1);

    return Intensity;
}
