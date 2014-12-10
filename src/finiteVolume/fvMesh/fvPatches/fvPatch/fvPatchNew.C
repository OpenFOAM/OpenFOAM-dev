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

#include "fvPatch.H"
#include "HashTable.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvPatch> Foam::fvPatch::New
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
{
    if (debug)
    {
        Info<< "fvPatch::New(const polyPatch&, const fvBoundaryMesh&) : "
            << "constructing fvPatch"
            << endl;
    }

    polyPatchConstructorTable::iterator cstrIter =
        polyPatchConstructorTablePtr_->find(patch.type());

    if (cstrIter == polyPatchConstructorTablePtr_->end())
    {
        FatalErrorIn("fvPatch::New(const polyPatch&, const fvBoundaryMesh&)")
            << "Unknown fvPatch type " << patch.type() << nl
            << "Valid fvPatch types are :"
            << polyPatchConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fvPatch>(cstrIter()(patch, bm));
}


// ************************************************************************* //
