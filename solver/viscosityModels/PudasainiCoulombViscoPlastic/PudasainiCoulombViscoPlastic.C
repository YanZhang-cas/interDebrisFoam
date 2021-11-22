/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
//howto compile new viscosity model: at folder incompressible (OpenFOAM/'username'-2.3.x/applications/transportModels)]$ wmake libso or wmake incompressible
/* rheology law with pressure dependend viscosity based on Pudasaini (2012): 'A general two-phase debris flow model'with 
   'delta' => internal friction angle, blending over to shear dependent viscosity due to 'my', 
   the higher 'my' the stronger the shear thinning. Using my = 1 will lead to nu*rhoS/sin(delta)/p = 1 at rest and declining in a way "between 
   logarhytmic and exponential" to about 0.1 at norm(D) reaching 10.

   Note about instability of incompressible Navier-Stokes equations with pressure-dependent viscosity in 
   'Some remarks on the Navier-Stokes equations with a pressure-dependent viscosity by Michael Renardy (1986) 
   published in 'Communications in partial differential equations':

   "Only as long as the eigenvalues of the symmetric part of the velocity gradient are less than 1/(2 nu') 
   [nu being the pressure dependent viscosity] can we expect to prove a local existence result.If this condition is violated, problems of 
   nonexistence and nonuniqueness occur in the constant coefficient problem and hence can be expected in the full Navier-Stokes system"

*/

#include "PudasainiCoulombViscoPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(PudasainiCoulombViscoPlastic, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        PudasainiCoulombViscoPlastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::PudasainiCoulombViscoPlastic::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    volScalarField limitedP = min(max(U_.db().lookupObject<volScalarField>("p_rgh"), dimensionedScalar ("PMIN", dimForce/dimLength/dimLength, VSMALL)), dimensionedScalar ("PMAX", dimForce/dimLength/dimLength, 1.e30));
    if (U_.db().foundObject<volScalarField>("p"))
	limitedP = min(max(U_.db().lookupObject<volScalarField>("p"), dimensionedScalar ("PMIN", dimForce/dimLength/dimLength, VSMALL)), dimensionedScalar ("PMAX", dimForce/dimLength/dimLength, 1.e30));

    volScalarField localAvLimP = max(fvc::average(fvc::interpolate(limitedP)),dimensionedScalar ("PMIN2", dimForce/dimLength/dimLength, VSMALL));
    // pressure oscillations possible due to pressure dependend viscosity, so do local averaging
/* Only as long as the eigenvalues of the symmetric part of the velocity gradient are less than 1/(2 nu') [nu being the pressure dependent viscosity] can we expect to prove a local existence result */
   // to get the max eigenvalues of the symmetric part of the velocity gradient, internal field and boundary field have to be adressed seperateley
    volVectorField eigSymmDU(eigenValues(symm(fvc::grad(U_))));
    scalarField maxEigSymmDU = cmptMax(eigSymmDU.internalField());
    // get from a vector field with vectors containing eigenvalues to a scalar field that only holds the largest eigenvalue:
    volScalarField magEigSymmDU = mag(eigSymmDU);
    magEigSymmDU.primitiveFieldRef() = maxEigSymmDU;

    /* limit at the boundary */
    // loop over all faces
	forAll(magEigSymmDU.boundaryFieldRef(), boundaryPatchID)
	    forAll(magEigSymmDU.boundaryFieldRef()[boundaryPatchID], faceI)
	    {
	       magEigSymmDU.boundaryFieldRef()[boundaryPatchID][faceI] = max(cmptMax(eigSymmDU.boundaryFieldRef()[boundaryPatchID][faceI]), 0.00000001);
	    }
    // now magEigSymmDU contains only a scalar. limit to values > 0 so that magEigSymmDU can be used in division:
    magEigSymmDU = max(magEigSymmDU, dimensionedScalar ("EigMIN", dimless/dimTime, 0.00000001));
	Info<< " calculate Coulomb-viscoplastic viscosity... " << endl;
/* pressure oscillations possible due to pressure dependend viscosity, so do local averaging*/
    volScalarField locAvNormD = calcLocAvNormD();
    volScalarField locAvNu = min( nu0_, nuMin_ + (sin(delta_)*localAvLimP/locAvNormD/rhoS_)*(1-exp(my_*(-1)*locAvNormD)));
    /* stability criteria asks for the derivative of nu with respect to pressure to be smaller than 1/(2*magEigSymmDU)/rhoS_. 
       Take a rough simplification that magEigSymmDU is pressure independent!*/
    volScalarField locNuLimit = localAvLimP/(2*magEigSymmDU)/rhoS_;
    Info<< " limit Coulomb-viscoplastic viscosity... " << endl;
    locAvNu = min(locAvNu, locNuLimit);
    locAvNu = max(locAvNu, nuMin_);
    return
    ( 
/* pressure dependent viscosity can cause local oscillation, average localy and limit to eigenvalues of the velocity gradient! */
        fvc::average(fvc::interpolate(locAvNu))
    );
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::PudasainiCoulombViscoPlastic::calcLocAvNormD() const
{
    tmp<volTensorField> DU(fvc::grad(U_)); //[1/sec]
    tmp<volTensorField> DUT(fvc::grad(U_)().T());
    tmp<volTensorField> DR(DU+DUT);
    tmp<volTensorField> D(0.5*DR); //[1/s]
    /* normD stores ||D|| 
	mag(D) = sqrt(D_ij D_ij ) = norm of strain rate tensor in OpenFOAM, where D_ij = 1/2 ((grad(U)+grad(U)^T) = symm(grad(U)) = 2*DR 
	the strain rates are a measure of how fast the three velocity components change in each of the three direction, three velocity components each vary in three directions => 3 x 3 strain rate tensor
	*/
    tmp<volScalarField> normD(max(mag(D), dimensionedScalar ("Dmin", dimless/dimTime, 0.0000000000001)));//[1/s]

    return
    ( 
/* pressure oscillations possible due to pressure dependend viscosity, so do local averaging! */
       fvc::average(fvc::interpolate(normD))
    );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::PudasainiCoulombViscoPlastic::PudasainiCoulombViscoPlastic
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    PudasainiCoulombViscoPlasticCoeffs_
    (
    	viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    delta_("delta", dimless, PudasainiCoulombViscoPlasticCoeffs_),
    rhoS_("rhoS", dimDensity, PudasainiCoulombViscoPlasticCoeffs_),
    nuMin_("nuMin", dimViscosity, PudasainiCoulombViscoPlasticCoeffs_),
    nu0_("nu0", dimViscosity, PudasainiCoulombViscoPlasticCoeffs_),
    my_("my", dimTime, PudasainiCoulombViscoPlasticCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    ),
    locAvNormD_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.mesh(), 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcLocAvNormD()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::PudasainiCoulombViscoPlastic::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    PudasainiCoulombViscoPlasticCoeffs_ = 
    		viscosityProperties.optionalSubDict(typeName + "Coeffs");

    PudasainiCoulombViscoPlasticCoeffs_.lookup("delta") >> delta_;
    PudasainiCoulombViscoPlasticCoeffs_.lookup("rhoS") >> rhoS_;
    PudasainiCoulombViscoPlasticCoeffs_.lookup("nuMin") >> nuMin_;
    PudasainiCoulombViscoPlasticCoeffs_.lookup("nu0") >> nu0_;
    PudasainiCoulombViscoPlasticCoeffs_.lookup("my") >> my_;

    return true;
}


// ************************************************************************* //
