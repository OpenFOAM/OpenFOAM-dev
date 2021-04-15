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

#include "psiThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiThermo, 0);
    defineRunTimeSelectionTable(psiThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiThermo> Foam::psiThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiThermo::~psiThermo()
{}


Foam::psiThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::psiThermo::correctRho(const Foam::volScalarField& deltaRho)
{}


Foam::tmp<Foam::volScalarField> Foam::psiThermo::implementation::rho() const
{
    return p()*psi_;
}


Foam::tmp<Foam::scalarField> Foam::psiThermo::implementation::rho
(
    const label patchi
) const
{
    return p().boundaryField()[patchi]*psi_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::psiThermo::implementation::rho0() const
{
    return p().oldTime()*psi_.oldTime();
}


const Foam::volScalarField& Foam::psiThermo::implementation::psi() const
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::psiThermo::implementation::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::psiThermo::implementation::mu
(
    const label patchi
) const
{
    return mu_.boundaryField()[patchi];
}


// ************************************************************************* //
