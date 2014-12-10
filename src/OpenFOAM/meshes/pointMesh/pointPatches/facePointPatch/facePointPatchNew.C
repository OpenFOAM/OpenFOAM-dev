/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "facePointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::facePointPatch> Foam::facePointPatch::New
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
{
    if (debug)
    {
        Info<< "facePointPatch::New(const polyPatch&, "
            << " const pointBoundaryMesh&) : "
            << "constructing facePointPatch"
            << endl;
    }

    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(patch.type());

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "facePointPatch::New(const polyPatch&, "
            "const pointBoundaryMesh&) : "
        )   << "Unknown facePointPatch type "
            << patch.type()
            << nl << nl
            << "Valid facePointPatch types are :" << endl
            << polyPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<facePointPatch>(cstrIter()(patch, bm));
}


// ************************************************************************* //
