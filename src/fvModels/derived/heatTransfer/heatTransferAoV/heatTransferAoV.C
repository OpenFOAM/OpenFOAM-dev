/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "heatTransferAoV.H"
#include "heatTransfer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferAoV::readCoeffs(const dictionary& dict)
{
    typeIOobject<volScalarField> AoVIO
    (
        "AoV",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (dict.found("AoV"))
    {
        AoV_ = dimensionedScalar("AoV", dimless/dimLength, dict);
        AoVPtr_.clear();
    }
    else if (AoVIO.headerOk())
    {
        AoV_ = dimensionedScalar("AoV", dimless/dimLength, NaN);
        AoVPtr_.set(new volScalarField(AoVIO, mesh_));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Area per unit volume (AoV) not found. A uniform AoV "
            << "value should be specified, or a non-uniform field should "
            << "exist at " << AoVIO.objectPath()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferAoV::heatTransferAoV
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    AoV_("AoV", dimless/dimLength, NaN),
    AoVPtr_(nullptr)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferAoV::~heatTransferAoV()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::heatTransferAoV::AoV() const
{
    if (!AoVPtr_.valid())
    {
        return volScalarField::New(typedName<heatTransfer>("AoV"), mesh_, AoV_);
    }
    else
    {
        return AoVPtr_();
    }
}


bool Foam::fv::heatTransferAoV::read(const dictionary& dict)
{
    readCoeffs(dict);

    return true;
}


// ************************************************************************* //
