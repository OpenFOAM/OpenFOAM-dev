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

#include "enginePiston.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enginePiston::enginePiston
(
    const fvMeshMover& meshMover,
    const word& pistonPatchName,
    const autoPtr<coordinateSystem>& pistonCS,
    const scalar minLayer,
    const scalar maxLayer
)
:
    meshMover_(refCast<const fvMeshMovers::engine>(meshMover)),
    patchIndex_(pistonPatchName, meshMover_.mesh().boundaryMesh()),
    csPtr_(pistonCS),
    minLayer_(minLayer),
    maxLayer_(maxLayer)
{}


Foam::enginePiston::enginePiston
(
    const fvMeshMover& meshMover,
    const dictionary& dict
)
:
    meshMover_(refCast<const fvMeshMovers::engine>(meshMover)),
    patchIndex_(dict.lookup("patch"), meshMover_.mesh().boundaryMesh()),
    csPtr_
    (
        coordinateSystem::New
        (
            meshMover_.mesh(),
            dict.subDict("coordinateSystem")
        )
    ),
    minLayer_(dict.lookup<scalar>("minLayer")),
    maxLayer_(dict.lookup<scalar>("maxLayer"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enginePiston::writeDict(Ostream& os) const
{
    os  << nl << token::BEGIN_BLOCK
        << "patch " << patchIndex_.name() << token::END_STATEMENT << nl
        << "minLayer " << minLayer_ << token::END_STATEMENT << nl
        << "maxLayer " << maxLayer_ << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //
