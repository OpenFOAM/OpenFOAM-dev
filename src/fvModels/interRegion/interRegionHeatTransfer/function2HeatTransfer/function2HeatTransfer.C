/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "function2HeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(function2HeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        function2HeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::function2HeatTransfer::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");
    UNbrName_ = coeffs().lookupOrDefault<word>("UNbr", "U");

    htcFunc_.reset(Function2<scalar>::New("htcFunc", coeffs()).ptr());
}


void Foam::fv::function2HeatTransfer::correctHtc() const
{
    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);

    const fvMesh& nbrMesh = mesh().time().lookupObject<fvMesh>(nbrRegionName());

    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);
    const scalarField UMagNbr(mag(UNbr));
    const scalarField UMagNbrMapped(interpolate(UMagNbr));

    htc_.primitiveFieldRef() = htcFunc_->value(mag(U()), UMagNbrMapped)*AoV_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::function2HeatTransfer::function2HeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    UName_(word::null),
    UNbrName_(word::null),
    htcFunc_(),
    AoV_
    (
        master()
      ? new volScalarField
        (
            IOobject
            (
                "AoV",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
      : nullptr
    )
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::function2HeatTransfer::~function2HeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::function2HeatTransfer::read(const dictionary& dict)
{
    if (interRegionHeatTransferModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
