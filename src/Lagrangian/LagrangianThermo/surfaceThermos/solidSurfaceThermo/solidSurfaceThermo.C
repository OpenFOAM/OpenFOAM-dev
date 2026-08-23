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

#include "solidSurfaceThermo.H"
#include "LagrangianSubMesh.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidSurfaceThermo, 0);
    defineRunTimeSelectionTable(solidSurfaceThermo, cloud);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::uniformDimensionedScalarField>
Foam::solidSurfaceThermo::implementation::p(const LagrangianSubMesh&) const
{
    return p_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidSurfaceThermo::implementation::implementation
(
    const basicLagrangianSubThermo& farThermo
)
:
    p_
    (
        IOobject
        (
            IOobject::groupName("p", farThermo.phaseName()),
            farThermo.mesh().time().name(),
            farThermo.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar
        (
            IOobject::groupName("p", farThermo.phaseName()),
            dimensions::pressure,
            NaN
        )
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidSurfaceThermo> Foam::solidSurfaceThermo::New
(
    const basicLagrangianSubThermo& farThermo
)
{
    return surfaceThermo::New<solidSurfaceThermo>(farThermo);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidSurfaceThermo::~solidSurfaceThermo()
{}


Foam::solidSurfaceThermo::implementation::~implementation()
{}


// ************************************************************************* //
