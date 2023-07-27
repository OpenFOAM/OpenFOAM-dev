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

#include "solidThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidThermo, 0);
    defineRunTimeSelectionTable(solidThermo, fvMesh);
}

const Foam::word Foam::solidThermo::derivedThermoName("heSolidThermo");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermo::implementation::implementation
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    p_
    (
        IOobject
        (
            phasePropertyName("p", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(phasePropertyName("p", phaseName), dimPressure, NaN)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermo::~solidThermo()
{}


Foam::solidThermo::implementation::~implementation()
{}


// ************************************************************************* //
