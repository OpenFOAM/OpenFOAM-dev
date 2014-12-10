/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const word& phaseName,
    const volScalarField& p,
    const volScalarField& T
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        p.mesh()
    ),
    name_(phaseName),
    p_(p),
    T_(T),
    thermo_(NULL),
    dgdt_
    (
        IOobject
        (
            IOobject::groupName("dgdt", phaseName),
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p.mesh(),
        dimensionedScalar("0", dimless/dimTime, 0)
    )
{
    {
        volScalarField Tp(IOobject::groupName("T", phaseName), T);
        Tp.write();
    }

    thermo_ = rhoThermo::New(p.mesh(), phaseName);
    thermo_->validate(phaseName, "e");

    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    notImplemented("phaseModel::clone() const");
    return autoPtr<phaseModel>(NULL);
}


void Foam::phaseModel::correct()
{
    thermo_->he() = thermo_->he(p_, T_);
    thermo_->correct();
}


// ************************************************************************* //
