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

#include "psiThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiThermo, 0);
    defineRunTimeSelectionTable(psiThermo, fvMesh);
}

const Foam::word Foam::psiThermo::derivedThermoName("hePsiThermo");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiThermo::implementation::implementation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName
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


Foam::tmp<Foam::volScalarField> Foam::psiThermo::renameRho()
{
    return rho();
}


Foam::tmp<Foam::volScalarField> Foam::psiThermo::implementation::rho() const
{
    return p()*psi();
}


Foam::tmp<Foam::scalarField> Foam::psiThermo::implementation::rho
(
    const label patchi
) const
{
    return p().boundaryField()[patchi]*psi().boundaryField()[patchi];
}


// ************************************************************************* //
