/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "heatTransferAv.H"
#include "heatTransfer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferAv::readCoeffs(const dictionary& dict)
{
    typeIOobject<volScalarField> AvIO
    (
        "Av",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    typeIOobject<volScalarField> AoVIO
    (
        "AoV",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (dict.found("Av"))
    {
        Av_ = dimensionedScalar("Av", dimless/dimLength, dict);
        AvPtr_.clear();
    }
    else if (dict.found("AoV"))
    {
        Av_ = dimensionedScalar("AoV", dimless/dimLength, dict);
        AvPtr_.clear();
    }
    else if (AvIO.headerOk())
    {
        Av_ = dimensionedScalar("Av", dimless/dimLength, NaN);
        AvPtr_.set(new volScalarField(AvIO, mesh_));
    }
    else if (AoVIO.headerOk())
    {
        Av_ = dimensionedScalar("AoV", dimless/dimLength, NaN);
        AvPtr_.set(new volScalarField(AoVIO, mesh_));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Area per unit volume (Av) not found. A uniform Av value "
            << "should be specified, or a non-uniform field should exist at "
            << AvIO.objectPath() << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferAv::heatTransferAv
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    Av_("Av", dimless/dimLength, NaN),
    AvPtr_(nullptr)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferAv::~heatTransferAv()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::heatTransferAv::Av() const
{
    if (!AvPtr_.valid())
    {
        return volScalarField::New(typedName<heatTransfer>("Av"), mesh_, Av_);
    }
    else
    {
        return AvPtr_();
    }
}


bool Foam::fv::heatTransferAv::read(const dictionary& dict)
{
    readCoeffs(dict);

    return true;
}


// ************************************************************************* //
