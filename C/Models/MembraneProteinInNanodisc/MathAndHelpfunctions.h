// Basic geometric operations
void SetRotationMatrix (double RM[3][3], double rz1, double rx, double rz2)
{
    RM[0][0] = cos(rz1) * cos(rz2) - sin(rz1) * cos(rx) * sin(rz2); RM[0][1] = -sin(rz1) * cos(rz2) - cos(rz1) * cos(rx) * sin(rz2); RM[0][2] =  sin(rx) * sin(rz2);
    RM[1][0] = cos(rz1) * sin(rz2) + sin(rz1) * cos(rx) * cos(rz2); RM[1][1] = -sin(rz1) * sin(rz2) + cos(rz1) * cos(rx) * cos(rz2); RM[1][2] = -sin(rx) * cos(rz2);
    RM[2][0] = sin(rz1) * sin(rx)                                 ; RM[2][1] =  cos(rz1) * sin(rx)                                 ; RM[2][2] =  cos(rx)           ;
}

void SetTranslationVector (double TV[3], double s, double theta, double z)
{
    TV[0] = s * cos(theta);
    TV[1] = s * sin(theta);
    TV[2] = z;
}


void Orient(double* x, double* y, double* z, double rotation[3][3], double translation[3])
{
    double xinit = *x;
    double yinit = *y;
    double zinit = *z;

    *x = rotation[0][0] * xinit + rotation[0][1] * yinit + rotation[0][2] * zinit + translation[0];
    *y = rotation[1][0] * xinit + rotation[1][1] * yinit + rotation[1][2] * zinit + translation[1];
    *z = rotation[2][0] * xinit + rotation[2][1] * yinit + rotation[2][2] * zinit + translation[2];
}

// Functions for handling complex arrays
double complex ** ComplexArray(int Dimension1, int Dimension2)
{
    int i;
	int j;
    double complex ** Array;

    Array = (double complex**) malloc(Dimension1 * sizeof(double complex*));

    for (i = 0; i < Dimension1; ++i) {
        Array[i] = (double complex*) malloc((Dimension2) * sizeof(double complex));

        for (j = 0; j < Dimension2; ++j) {
            Array[i][j] = 0.0;
        }
    }

    return Array;
}

void FreeComplexArray(double complex **Alpha, int Dimension1, int Dimension2)
{
    int i;

    for (i = 0; i < Dimension1; ++i) {
        free(Alpha[i]);
    }

    free(Alpha);
}

//// Function for copying residue
//void CopyResidue(struct Residue * Original, struct Residue * Copy)
//{
//      Copy->xVolume = Original->xVolume;
//      Copy->yVolume = Original->yVolume;
//      Copy->zVolume = Original->zVolume;
//      
//      Copy->xXRayScattering = Original->xXRayScattering;
//      Copy->yXRayScattering = Original->yXRayScattering;
//      Copy->zXRayScattering = Original->zXRayScattering;
//      
//      Copy->xNeutronScattering = Original->xNeutronScattering;
//      Copy->yNeutronScattering = Original->yNeutronScattering;
//      Copy->zNeutronScattering = Original->zNeutronScattering;
//    
//      Copy->XRayScatteringLength = Original->XRayScatteringLength;
//      Copy->NeutronScatteringLength = Original->NeutronScatteringLength;
//      Copy->Volume = Original->Volume;

//     Copy->Name[0] = Original->Name[0];
//     Copy->Name[1] = Original->Name[1];
//     Copy->Name[2] = Original->Name[2];
//     Copy->ResidueID = Original->ResidueID;
//}


// changed by MS
void CopyResidue(struct Residue * Original, struct Residue * Copy)
{
	Copy->xVolume = Original->xVolume;
	Copy->yVolume = Original->yVolume;
	Copy->zVolume = Original->zVolume;

	Copy->xXRayScattering = Original->xXRayScattering;
	Copy->yXRayScattering = Original->yXRayScattering;
	Copy->zXRayScattering = Original->zXRayScattering;

	Copy->xNeutronScattering = Original->xNeutronScattering;
	Copy->yNeutronScattering = Original->yNeutronScattering;
	Copy->zNeutronScattering = Original->zNeutronScattering;

	Copy->XRayScatteringLength    = Original->XRayScatteringLength;
	Copy->NeutronScatteringLength = Original->NeutronScatteringLength;

	Copy->Volume    = Original->Volume;
	Copy->Weight    = Original->Weight ; // added by MS

	Copy->Name[0]   = Original->Name[0];
	Copy->Name[1]   = Original->Name[1];
	Copy->Name[2]   = Original->Name[2];
	Copy->Name[3]   = Original->Name[3];

	Copy->AtomName[0]   = Original->AtomName[0]; // added by MS
	Copy->AtomName[1]   = Original->AtomName[1]; // added by MS
	Copy->AtomName[2]   = Original->AtomName[2]; // added by MS

	Copy->ResidueID = Original->ResidueID;
}


// Functions for handling complex numbers
double Sinc(double x) {
    double Result;

    if (x == 0.0) {
        Result = 1.0;
    } else {
        Result = sin(x) / x;
	}

    return Result;
}

int Sign(double x) {

    if (x < 0.0) {
        return -1;
    } else if (x > 0.0) {
        return 1;
    } else {
        return 0;
	}
}

double complex pol(double r, double phi)
{
    if (phi == 0.0) {
        return r + I * 0.0;
    } else {
        return r * (cos(phi) + I * sin(phi));
	}
}

// Expand analytical formfactors on Alpha - from Søren Kynde (2011) - Martin Cramer Pedersen (2017)
void ExpansionOnSphericalHarmonics(double complex ** Alpha, double complex ** DummyArray)
{
	int i;
	int j;
    int l;
	int m;

    int ntheta;
	int nphi;

    double theta;
	double phi;
    double thetastep = M_PI / Ntheta;
	double phistep = 2.0 * M_PI / Nphi;

    double complex fm[Ntheta];
    double complex phase[Nh + 1][Nphi];
    double sinth[Ntheta];
    double Weights[Ntheta];
    double Legendre[Ntheta][Nh + 1];

	// Zero arrays
    for (i = 0; i < Ntheta; ++i) {
        Weights[i] = 0.0;
		fm[i] = 0.0;
    }

	// Weights for theta integral
    for (j = 0; j < Ntheta; ++j) { 
        theta = (j + 0.5) * thetastep;

        for (l = 0; l < Ntheta / 2; ++l) {
            Weights[j] += 2.0 / (Ntheta / 2) * 1.0 / (2 * l + 1) * sin((2 * l + 1) * theta);
        }
    }

	// Skip every second spherical harmonic due to the symmetry of the disc
    for (m = 0; m < Nh + 1; m += skip) {

        for (ntheta = 0; ntheta < Ntheta; ++ntheta) {
            fm[ntheta] = 0.0;

            for (nphi = 0; nphi < Nphi; ++nphi) {
                phi = phistep * (nphi + 0.5);

                if (ntheta == 0) {
                    phase[m][nphi] = pol(phistep, -m * phi);
				}

				// fm(theta) = int_0^2pi DummyArray(theta, phi) exp(-m * phi) dphi
                fm[ntheta] += DummyArray[ntheta][nphi] * phase[m][nphi];
            }
        }

        for (l = m; l < Nh + 1; l += skip) {

            for (ntheta = 0; ntheta < Ntheta; ++ntheta) {
                theta = (ntheta + 0.5) * thetastep;

                if (l == m) {
                    gsl_sf_legendre_sphPlm_array(Nh, m, cos(theta), &Legendre[ntheta][l]);
                    sinth[ntheta] = sin(theta);
                }

				// Alpha_lm = int_0^pi P_lm sin(theta) fm(theta) dtheta
                Alpha[l][m] += 1.0 / sqrt(4.0 * M_PI) * Legendre[ntheta][l] * Weights[ntheta] * sinth[ntheta] * fm[ntheta];
            }
        }
    }
}


bool ResidueIsInLayer(double z, double Height)
{
	return (fabs(z) < Height / 2.0);
}


// Scattering from residues - Søren Kynde (2012) - Martin Pedersen (2017)
void AddScatteringFromResidue(double complex ** Beta, double q, struct Residue CurrentResidue, double Contrast, double SLDOfHeads, double SLDOfAlkyls, double SLDOfMethyls, double SLDOfBelt, double SLDOfSolvent, double HeightOfLipids, double HeightOfCore, double HeightOfMethyl, double HeightOfBelt, double CorrectionToVolumeOfMP, double MajorSemiaxisOfCore, double MinorSemiaxisOfCore)
{
    int l;
	int m;
    double x;
    double y;
    double z;
    double Legendre[Nh + 1];
    double Bessel[Nh + 1];
    double SLDOfBackground;
    double ScatteringLengthOfResidue;
	double ScatteringLengthOfSolvent;
	double DifferenceInScatteringLength;

	// Determine background
	SLDOfBackground = SLDOfSolvent;

    if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfMethyl)) {
        SLDOfBackground = SLDOfMethyls;
    } else if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfCore)) {
        SLDOfBackground = SLDOfAlkyls;
    } else if (ResidueIsInLayer(CurrentResidue.zVolume, HeightOfLipids)) {
        SLDOfBackground = SLDOfHeads;
	}

	// Calculation location of bead
    if (Contrast < 0.0 || Contrast > 100.0) {
       ScatteringLengthOfResidue    = CurrentResidue.XRayScatteringLength;
       ScatteringLengthOfSolvent    = CurrentResidue.Volume * SLDOfBackground * CorrectionToVolumeOfMP;
       DifferenceInScatteringLength = ScatteringLengthOfResidue - ScatteringLengthOfSolvent;

       x = (CurrentResidue.xXRayScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
       y = (CurrentResidue.yXRayScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
       z = (CurrentResidue.zXRayScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
    } else {
       ScatteringLengthOfResidue    = CurrentResidue.NeutronScatteringLength;
       ScatteringLengthOfSolvent    = CurrentResidue.Volume * SLDOfBackground * CorrectionToVolumeOfMP;
       DifferenceInScatteringLength = ScatteringLengthOfResidue - ScatteringLengthOfSolvent;

       x = (CurrentResidue.xNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.xVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
       y = (CurrentResidue.yNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.yVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
       z = (CurrentResidue.zNeutronScattering * ScatteringLengthOfResidue - CurrentResidue.zVolume * ScatteringLengthOfSolvent) / DifferenceInScatteringLength;
    }

	// Change to spherical coordinates
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double theta = acos(z / r);
    double phi = acos(x / (r * sin(theta))) * Sign(y);

	// Spherical bessel functions
    gsl_sf_bessel_jl_array(Nh, q * r, Bessel);

    for (m = 0; m < Nh + 1; ++m) {
		// Legendre polynomials
        gsl_sf_legendre_sphPlm_array(Nh, m, cos(theta), &Legendre[m]);
                                                                    
        for (l = m; l < Nh + 1; ++l) {
            Beta[l][m] += sqrt(4.0 * M_PI) * cpow(I, l) * DifferenceInScatteringLength * Bessel[l] * Legendre[l] * pol(1.0, -m * phi);
        }
    }
}

// The form factor of a disc in matrix form - Søren Kynde (2012) - Martin Pedersen (2017)
void Discflat(double complex ** DummyArray, double MajorSemiaxis, double MinorSemiaxis, double Height, double SLD, double r0, double theta0, double phi0, double q)
{
    int ntheta;
	int nphi;
    double thetastep = M_PI / Ntheta;
	double phistep = 2.0 * M_PI / Nphi;

    double cosr0q = 0.0;
    double sinc;

    double r;
	double theta;
	double phi;

    double Volume = M_PI * MajorSemiaxis * MinorSemiaxis * Height;

    for (ntheta = 0; ntheta < Ntheta; ++ntheta) {
        theta = (ntheta + 0.5) * thetastep;
        sinc = Sinc(Height / 2.0 * q * cos(theta));

        for (nphi = 0; nphi < Nphi; ++nphi) {
            phi = (nphi + 0.5) * phistep;

            if ((r0 != 0.0) || (theta0 != 0.0) || (phi0 != 0.0)) {
                cosr0q = sin(theta0) * sin(theta) * cos(phi0 - phi) + cos(theta0) * cos(theta);
			}

            r = sin(theta) * sqrt(pow(MajorSemiaxis * cos(phi), 2) + pow(MinorSemiaxis * sin(phi), 2));

            if (q * r == 0.0) {
                DummyArray[ntheta][nphi] += Volume * SLD * pol(1.0 * sinc, 0.0);
            } else {
                DummyArray[ntheta][nphi] += Volume * SLD * pol(2.0 * j1(q * r) / (q * r) * sinc, -r0 * q * cosr0q);
			}
        }
    }
}

// Collection of amplitudes
void NanodiscPDBModel(complex double ** Alpha, double q, double SLDOfHeads, double SLDOfAlkyls, double SLDOfMethyls, double SLDOfBelt, double SLDOfSolvent, double MajorSemiaxisOfCore, double MinorSemiaxisOfCore, double ThicknessOfBelt, double HeightOfLipids, double HeightOfCore, double HeightOfMethyl, double HeightOfBelt)
{
    double complex ** DummyArray = ComplexArray(Ntheta, Nphi);

	// Calculate formfactor of Nanodisc consisting of concentric eliptical cylinders
    Discflat(DummyArray, MajorSemiaxisOfCore                  , MinorSemiaxisOfCore                  , HeightOfLipids, SLDOfHeads - SLDOfSolvent , 0.0, 0.0, 0.0, q);
    Discflat(DummyArray, MajorSemiaxisOfCore                  , MinorSemiaxisOfCore                  , HeightOfCore  , SLDOfAlkyls - SLDOfHeads  , 0.0, 0.0, 0.0, q);
    Discflat(DummyArray, MajorSemiaxisOfCore                  , MinorSemiaxisOfCore                  , HeightOfMethyl, SLDOfMethyls - SLDOfAlkyls, 0.0, 0.0, 0.0, q);
    Discflat(DummyArray, MajorSemiaxisOfCore + ThicknessOfBelt, MinorSemiaxisOfCore + ThicknessOfBelt, HeightOfBelt  , SLDOfBelt - SLDOfSolvent  , 0.0, 0.0, 0.0, q);
    Discflat(DummyArray, MajorSemiaxisOfCore                  , MinorSemiaxisOfCore                  , HeightOfBelt  , SLDOfSolvent - SLDOfBelt  , 0.0, 0.0, 0.0, q);
   
    // Expand formfactor amplitud in spherical harmonics with coefficients Alpha
    ExpansionOnSphericalHarmonics(Alpha, DummyArray);
    FreeComplexArray(DummyArray, Ntheta, Nphi);
}
