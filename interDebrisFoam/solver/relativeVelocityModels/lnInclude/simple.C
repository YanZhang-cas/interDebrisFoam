/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "simple.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(simple, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, simple, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simple::simple
(
    const dictionary& dict,
    const incompressibleThreePhaseMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    a_("a", dimless, dict),
    V0_("V0", dimVelocity, dict),
    d_("d", dimLength, dict),
    g_("g",dimAcceleration,dict),
    residualAlpha_("residualAlpha", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::simple::~simple()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::simple::correct()
{
    //V0_ = d_*d_*pow(max(alphac_, scalar(0)), 4.0)*(rhod()-rhoc())*9.81/(18.0*rhoc()*mixture().nuModel3_->nu());
    //U23_ = V0_*pow(max(alphad_, scalar(0)), a_);
    //volScalarField muc = max(alphac_,scalar(0.0))*rhoc_*mixture().nuModel3().nu();
    //volScalarField mud = max(alphad_,scalar(0.0))*rhod_*mixture().nuModel2().nu();
    //dimensionedScalar small("0.001", dimViscosity*dimDensity, 0.00001);
    
    //U23_ = d_*d_*max(alphad_,scalar(0.0))*(rhod_-rhoc_)*g_/(18.0*(muc + mud + small));
    U23_=V0_*pow(max(1.0-alphad_,scalar(0.0)), a_);
}


// ************************************************************************* //
