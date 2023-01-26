/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "function2.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferCoefficientModels
{
    defineTypeNameAndDebug(function2, 0);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, function2, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::function2::readCoeffs
(
    const dictionary& dict
)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");
    UNbrName_ = dict.lookupOrDefault<word>("UNbr", "U");

    htcFunc_.reset(Function2<scalar>::New("htcFunc", dict).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::function2::function2
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    heatTransferCoefficientModel(typeName, dict, model),
    model_(model),
    UName_(word::null),
    UNbrName_(word::null),
    htcFunc_(),
    htc_
    (
        IOobject
        (
            typedName("htc"),
            model.mesh().time().name(),
            model.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        model.mesh(),
        dimensionedScalar(dimPower/dimTemperature/dimArea, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::function2::~function2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::function2::correct()
{
    const scalarField UMag
    (
        mag
        (
            model_.mesh()
           .lookupObject<volVectorField>(UName_)
           .primitiveField()
        )
    );
    const scalarField UMagNbr
    (
        model_.interpolate
        (
            mag
            (
                model_.nbrMesh()
               .lookupObject<volVectorField>(UNbrName_)
               .primitiveField()
            )()
        )
    );

    htc_.primitiveFieldRef() = htcFunc_->value(UMag, UMagNbr);
    htc_.correctBoundaryConditions();
}


bool Foam::fv::heatTransferCoefficientModels::function2::read
(
    const dictionary& dict
)
{
    if (heatTransferCoefficientModel::read(dict))
    {
        readCoeffs(dict);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
