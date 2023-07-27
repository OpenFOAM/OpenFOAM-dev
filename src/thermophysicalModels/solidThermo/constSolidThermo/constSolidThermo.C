/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "constSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(constSolidThermo, 0);
    addToRunTimeSelectionTable(basicThermo, constSolidThermo, fvMesh);
    addToRunTimeSelectionTable(solidThermo, constSolidThermo, fvMesh);
}


// * * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * //

Foam::constSolidThermo::constSolidThermo
(
    const fvMesh& mesh,
    const bool readKappa,
    const word& phaseName
)
:
    PhysicalPropertiesThermo<solidThermo::composite>(mesh, phaseName),
    Cv_(readProperty<scalar>("Cv", dimEnergy/dimMass/dimTemperature)),
    e_
    (
        IOobject
        (
            phasePropertyName("e"),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Cv_*T_,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    rho_ = readProperty<scalar>("rho", rho_.dimensions());

    if (readKappa)
    {
        kappa_ = readProperty<scalar>("kappa", kappa_.dimensions());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSolidThermo::constSolidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermo(mesh, true, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constSolidThermo::~constSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::constSolidThermo::Cpv() const
{
    return Cv_;
}


Foam::volScalarField& Foam::constSolidThermo::he()
{
    return e_;
}


const Foam::volScalarField& Foam::constSolidThermo::he() const
{
    return e_;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return scalarField(Cv_, cells)*T;
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return Cv_*T;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi]*T;
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hs() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::ha() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::hc() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::THE
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::THE
(
    const scalarField& he,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


const Foam::volScalarField& Foam::constSolidThermo::Cp() const
{
    return Cv_;
}


const Foam::volScalarField& Foam::constSolidThermo::Cv() const
{
    return Cv_;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return Cv_.boundaryField()[patchi];
}


const Foam::volVectorField& Foam::constSolidThermo::Kappa() const
{
    NotImplemented;
    return volVectorField::null();
}


void Foam::constSolidThermo::correct()
{
    T_ = e_/Cv_;
}


// ************************************************************************* //
