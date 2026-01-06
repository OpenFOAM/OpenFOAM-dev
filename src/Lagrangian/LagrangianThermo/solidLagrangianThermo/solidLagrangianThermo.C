/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "solidLagrangianThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidLagrangianThermo, 0);
    defineRunTimeSelectionTable(solidLagrangianThermo, LagrangianMesh);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::uniformDimensionedScalarField>
Foam::solidLagrangianThermo::implementation::p(const LagrangianSubMesh&) const
{
    return p_;
}


Foam::tmp<Foam::uniformDimensionedScalarField>
Foam::solidLagrangianThermo::implementation::p
(
    const LagrangianInjection&,
    const LagrangianSubMesh&
) const
{
    return p_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidLagrangianThermo::implementation::implementation
(
    const dictionary& dict,
    const LagrangianMesh& mesh,
    const word& phaseName
)
:
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(IOobject::groupName("p", phaseName), dimPressure, NaN)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidLagrangianThermo>
Foam::solidLagrangianThermo::New
(
    const LagrangianMesh& mesh,
    const word& phaseName
)
{
    return basicLagrangianThermo::New<solidLagrangianThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidLagrangianThermo::~solidLagrangianThermo()
{}


Foam::solidLagrangianThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solidLagrangianThermo::initialise()
{
    correct(mesh().subAll());
}


void Foam::solidLagrangianThermo::correctPressure(const LagrangianSubMesh&)
{}


// ************************************************************************* //
