/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "compressibleVoFphase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleVoFphase::compressibleVoFphase
(
    const word& name,
    const fvMesh& mesh,
    const volScalarField& T
)
:
    VoFphase(name, mesh),
    thermo_(nullptr),
    Alpha_
    (
        IOobject
        (
            IOobject::groupName("Alpha", name),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    ),
    vDot_
    (
        IOobject
        (
            IOobject::groupName("vDot", name),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, 0)
    )
{
    {
        volScalarField Tp(IOobject::groupName("T", name), T);
        Tp.write();
    }

    thermo_ = rhoFluidThermo::New(mesh, name);
    thermo_->validate(name, "e");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::VoFphase> Foam::compressibleVoFphase::clone() const
{
    NotImplemented;
    return autoPtr<VoFphase>(nullptr);
}


void Foam::compressibleVoFphase::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    thermo_->he() = thermo_->he(p, T);
    thermo_->correct();
}


// ************************************************************************* //
