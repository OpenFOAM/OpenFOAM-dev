/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
        option,
        function2HeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Function2<Foam::scalar>&
Foam::fv::function2HeatTransfer::htcFunc() const
{
    if (!htcFunc_.valid())
    {
        htcFunc_ = Function2<scalar>::New("htcFunc", coeffs_);
    }

    return htcFunc_();
}


const Foam::volScalarField& Foam::fv::function2HeatTransfer::AoV() const
{
    if (!AoV_.valid())
    {
        AoV_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "AoV",
                    startTimeName_,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    return AoV_();
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
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    UNbrName_(coeffs_.lookupOrDefault<word>("UNbr", "U")),
    htcFunc_(),
    AoV_(),
    startTimeName_(mesh.time().timeName())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::function2HeatTransfer::~function2HeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::function2HeatTransfer::calculateHtc() const
{
    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);

    const scalarField UMagNbr(mag(UNbr));
    const scalarField UMagNbrMapped(interpolate(UMagNbr));
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    htc_.primitiveFieldRef() = htcFunc().value(mag(U()), UMagNbrMapped)*AoV();
}


bool Foam::fv::function2HeatTransfer::read(const dictionary& dict)
{
    if (interRegionHeatTransferModel::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
