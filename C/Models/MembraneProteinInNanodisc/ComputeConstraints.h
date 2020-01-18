void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    // Declare dummy variables needed in function
    double AreaOfDisc;
	int i = 0;
    int j;

    // Variables describing the protein orientation
    double VA;
    double VM;
    double VH;

    struct Residue CurrentResidue;
    double TranslationVector[3];
    double RotationMatrix[3][3];
    SetTranslationVector(TranslationVector, Parameters[24], Parameters[25], Parameters[26]);
    SetRotationMatrix(RotationMatrix, Parameters[27], Parameters[28], Parameters[29]);

    // Variables describing the lipid core
    double HeightOfCore;
    double HeightOfBelt;
    double HeightOfLipids;
    double HeightOfMethyl;

    // Variables describing scattering lengths
    double ScatteringLengthDensityOfCaps;
    double ScatteringLengthDensityOfCore;
    double ScatteringLengthDensityOfBelt;
    double ScatteringLengthDensityOfWater;
    double ScatteringLengthDensityOfMethyl;

    double ScatteringLengthOfBelt;
    double ScatteringLengthOfCore;
    double ScatteringLengthOfWater;
    double ScatteringLengthOfCaps;
    double ScatteringLengthOfMethyl;

    // Variables describing the volumes
    double VolumeOfWater;
    double VolumeOfHead;
    double VolumeOfHeads;
    double VolumeOfTail;
    double VolumeOfTails;
    double VolumeOfMethyl;
    double VolumeOfMethyls;
    double VolumeOfBelt;

    double CorrectionToVolumeOfLipid;
    double CorrectionToVolumeOfWater;
    double CorrectionToVolumeOfBelt;
    double CorrectionToVolumeOfMP;

    // Variables describing water
    double NumberOfWaterAtBelt;
    double NumberOfWaterAtHead;

    // Variables describing the properties of the disc
    double RatioBetweenAxis;
    double ThicknessOfBelt;

    double SemiMajorAxisOfCore;
    double SemiMinorAxisOfCore;

    // Variables describing the lipids
    double NumberOfLipids;
    double DisplacedAlkyl0;
    double DisplacedMethyl0;
    double DisplacedAlkyl  = 10.0;
    double DisplacedMethyl = 10.0;
    double DisplacedHeads  = 10.0;

    /// Get parameters from Parameters
    CorrectionToVolumeOfBelt  = Parameters[9];
    CorrectionToVolumeOfLipid = Parameters[10];
    CorrectionToVolumeOfWater = Parameters[19];
    CorrectionToVolumeOfMP    = Parameters[11];

    // Hydration numbers
    NumberOfWaterAtHead = fabs(Parameters[5]);
    NumberOfWaterAtBelt = Parameters[6];

    // Get parameters from variable VolumesOfMolecules
    VolumeOfWater  =       VolumesOfMolecules[0];
    VolumeOfHead   =       VolumesOfMolecules[1] * CorrectionToVolumeOfLipid + NumberOfWaterAtHead * VolumesOfMolecules[0] * CorrectionToVolumeOfWater;
    VolumeOfTail   =       VolumesOfMolecules[2] * CorrectionToVolumeOfLipid;
    VolumeOfMethyl =       VolumesOfMolecules[3] * CorrectionToVolumeOfLipid;
    VolumeOfBelt   = 2.0 * VolumesOfMolecules[4] * CorrectionToVolumeOfBelt + 2.0 * NumberOfWaterAtBelt * VolumesOfMolecules[0];

    // Obtain scattering lengths
    ScatteringLengthOfWater  =       ScatteringLengths[0];
    ScatteringLengthOfCaps   =       ScatteringLengths[1] + NumberOfWaterAtHead * ScatteringLengths[0];
    ScatteringLengthOfCore   =       ScatteringLengths[2];
    ScatteringLengthOfMethyl =       ScatteringLengths[3];
    ScatteringLengthOfBelt   = 2.0 * ScatteringLengths[4] + 2.0 * NumberOfWaterAtBelt * ScatteringLengths[0];

    /// Derive scattering lengths
    ScatteringLengthDensityOfWater  = ScatteringLengthOfWater  / VolumeOfWater;
    ScatteringLengthDensityOfCaps   = ScatteringLengthOfCaps   / VolumeOfHead;
    ScatteringLengthDensityOfCore   = ScatteringLengthOfCore   / VolumeOfTail;
    ScatteringLengthDensityOfMethyl = ScatteringLengthOfMethyl / VolumeOfMethyl;
    ScatteringLengthDensityOfBelt   = ScatteringLengthOfBelt   / VolumeOfBelt;
        
    // Assign values of the constraints
    Constraints[0] = ScatteringLengthDensityOfCaps;
    Constraints[1] = ScatteringLengthDensityOfCore;
    Constraints[2] = ScatteringLengthDensityOfMethyl;
    Constraints[3] = ScatteringLengthDensityOfBelt;
    Constraints[4] = ScatteringLengthDensityOfWater;

	// Geometric parameters
    RatioBetweenAxis = Parameters[0];
    HeightOfBelt     = Parameters[2];
    double AreaPerHeadgroup = Parameters[1];

     HeightOfCore      = 2.0 * (VolumeOfTail + VolumeOfMethyl) / AreaPerHeadgroup;

    // Number of lipids
    NumberOfLipids = Parameters[3];

    do{
        DisplacedAlkyl0  = DisplacedAlkyl;
        DisplacedMethyl0 = DisplacedMethyl;

        // Volume of the lipid layers including the parts displaced by the membrane protein
        VolumeOfHeads   = VolumeOfHead   * (NumberOfLipids + DisplacedHeads);
        VolumeOfTails   = VolumeOfTail   * (NumberOfLipids + DisplacedAlkyl); 
        VolumeOfMethyls = VolumeOfMethyl * (NumberOfLipids + DisplacedMethyl);

        HeightOfMethyl      = HeightOfCore * VolumeOfMethyls / (VolumeOfTails + VolumeOfMethyls);
        SemiMinorAxisOfCore = sqrt(VolumeOfMethyls / (HeightOfMethyl * pi * RatioBetweenAxis));
        SemiMajorAxisOfCore = SemiMinorAxisOfCore * RatioBetweenAxis;

        AreaOfDisc = pi * SemiMinorAxisOfCore * SemiMajorAxisOfCore;
        HeightOfLipids = VolumeOfHeads / AreaOfDisc + HeightOfCore;

        VA = 0.0;
		VM = 0.0;
		VH = 0.0;

        for (j = 0; j < ProteinStructure.NumberOfResidues; ++j) {
            CopyResidue(&ProteinStructure.Residues[j], &CurrentResidue);
            Orient(&CurrentResidue.xVolume, &CurrentResidue.yVolume, &CurrentResidue.zVolume, RotationMatrix, TranslationVector);

            if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfMethyl)) {
                VM += CurrentResidue.Volume * CorrectionToVolumeOfMP;
            } else if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfCore)) {
                VA += CurrentResidue.Volume * CorrectionToVolumeOfMP;
            } else if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfLipids)) {
                VH += CurrentResidue.Volume * CorrectionToVolumeOfMP;
			}
        }

        DisplacedHeads  = VH / VolumeOfHead;
        DisplacedAlkyl  = VA / VolumeOfTail;
        DisplacedMethyl = VM / VolumeOfMethyl;

        ++i;

    } while (fabs(DisplacedMethyl - DisplacedMethyl0) / DisplacedMethyl + fabs(DisplacedAlkyl - DisplacedAlkyl0) / DisplacedAlkyl > 0.1 && i < 20);

    // Volume of the lipid layers including the parts displaced by the membrane protein
    VolumeOfHeads   = VolumeOfHead   * (NumberOfLipids + DisplacedHeads);
    VolumeOfTails   = VolumeOfTail   * (NumberOfLipids + DisplacedAlkyl); 
    VolumeOfMethyls = VolumeOfMethyl * (NumberOfLipids + DisplacedMethyl);

    HeightOfMethyl      = HeightOfCore * VolumeOfMethyls / (VolumeOfTails + VolumeOfMethyls);
    SemiMinorAxisOfCore = sqrt(VolumeOfMethyls / (HeightOfMethyl * pi * RatioBetweenAxis));
    SemiMajorAxisOfCore = SemiMinorAxisOfCore * RatioBetweenAxis;
    
    AreaOfDisc      = pi * SemiMinorAxisOfCore * SemiMajorAxisOfCore;
    HeightOfLipids  = VolumeOfHeads / AreaOfDisc + HeightOfCore;
    ThicknessOfBelt = - (SemiMajorAxisOfCore + SemiMinorAxisOfCore) / 2.0 + sqrt(pow(SemiMajorAxisOfCore + SemiMinorAxisOfCore, 2) / 4.0 + VolumeOfBelt / (pi * HeightOfBelt));

  

    // Assign values of the constraints
    Constraints[6] = HeightOfLipids;
    Constraints[7] = HeightOfCore;
    Constraints[8] = HeightOfMethyl;
    Constraints[21] = Concentration;
    Constraints[23] = SemiMajorAxisOfCore;
    Constraints[24] = SemiMinorAxisOfCore;
    Constraints[27] = ThicknessOfBelt;
    Constraints[9] = AreaPerHeadgroup;
}
