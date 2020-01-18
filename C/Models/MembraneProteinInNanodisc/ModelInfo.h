/************************************************
 *
 * Model name: Membrane protein in nanodisc
 * Author    : Martin Cramer Pedersen (based on code by Søren Kynde)
 * Email     : mcpe@nbi.ku.dk
 *
 ***********************************************/

// List of headers to be included on compilation
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#define Nh 19
#define skip 2
#define Nphi ((Nh+1)*2)
#define Ntheta ((Nh+1)*2)

#include "MathAndHelpfunctions.h"
#include "ComputeConstraints.h"
#include "Model.h"
#include "OutputData.h"
