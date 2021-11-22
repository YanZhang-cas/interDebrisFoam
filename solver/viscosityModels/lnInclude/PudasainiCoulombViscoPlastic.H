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

Class
    Foam::viscosityModels::PudasainiCoulombViscoPlastic

Description
     Visco-Palstic non-Newtonian viscosity model.

SourceFiles
    PudasainiCoulombViscoPlastic.C

\*---------------------------------------------------------------------------*/

#ifndef PudasainiCoulombViscoPlastic_H
#define PudasainiCoulombViscoPlastic_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "fvc.H"
#include "tensor.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class PudasainiCoulombViscoPlastic Declaration
\*---------------------------------------------------------------------------*/

class PudasainiCoulombViscoPlastic
:
    public viscosityModel
{
    // Private data

        dictionary PudasainiCoulombViscoPlasticCoeffs_;

        dimensionedScalar delta_;
        dimensionedScalar rhoS_;
        dimensionedScalar nuMin_;
        dimensionedScalar nu0_;
        dimensionedScalar my_;

        volScalarField nu_;
        volScalarField locAvNormD_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        tmp<volScalarField> calcLocAvNormD() const;

public:

    //- Runtime type information
    TypeName("PudasainiCoulombViscoPlastic");


    // Constructors

        //- Construct from components
        PudasainiCoulombViscoPlastic
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~PudasainiCoulombViscoPlastic()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }
	
        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
	
        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
	    locAvNormD_ = calcLocAvNormD();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
