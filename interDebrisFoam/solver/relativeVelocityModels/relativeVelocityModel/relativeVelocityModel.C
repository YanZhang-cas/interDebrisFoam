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

#include "relativeVelocityModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions   * * * * * * * * * * * //

Foam::wordList Foam::relativeVelocityModel::U23PatchFieldTypes() const
{
    const volVectorField& U = mixture_.U();

    wordList U23Types
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    forAll(U.boundaryField(), i)
    {
        if
        (
            isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
         || isA<slipFvPatchVectorField>(U.boundaryField()[i])
         || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
        )
        {
            U23Types[i] = fixedValueFvPatchVectorField::typeName;
        }
    }

    return U23Types;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleThreePhaseMixture& mixture
)
:
    mixture_(mixture),
    alphac_(mixture.alpha2()),
    alphad_(mixture.alpha3()),
    rhoc_(mixture.rho2()),
    rhod_(mixture.rho3()),

    U23_
    (
        IOobject
        (
            "U23",
            alphac_.time().timeName(),
            alphac_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alphac_.mesh(),
        dimensionedVector("U23", dimVelocity, Zero),
        U23PatchFieldTypes()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleThreePhaseMixture& mixture
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown time scale model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid time scale model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            cstrIter()
            (
                dict.optionalSubDict(modelType + "Coeffs"),
                mixture
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> Foam::relativeVelocityModel::rho() const
{
    return alphad_*rhod_ + alphac_*rhoc_;
}

//tmp<dimensionedScalar> Foam::relativeVelocityModel::rhoc()
//{
//    return rhoc_;
//}

//tmp<dimensionedScalar> Foam::relativeVelocityModel::rhod()
//{
//    return rhod_;
//}


tmp<volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    dimensionedScalar small("0", dimDensity, scalar(0.00001));

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            "tauDm",
            max(alphad_,scalar(0.0))*max(alphac_,scalar(0.0))*rhoc_*rhod_*sqr(U23_)/(rho()+small)
        )
    );
}


// ************************************************************************* //
