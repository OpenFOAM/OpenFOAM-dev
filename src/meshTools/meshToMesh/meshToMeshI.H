/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "meshToMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline const Foam::List<Foam::labelPair>&
Foam::meshToMesh::srcToTgtPatchIDs() const
{
    return srcToTgtPatchIDs_;
}


inline const Foam::PtrList<Foam::patchToPatch>&
Foam::meshToMesh::srcToTgtPatchToPatches() const
{
    return srcToTgtPatchToPatches_;
}


inline const Foam::labelList& Foam::meshToMesh::tgtCuttingPatchIDs() const
{
    return tgtCuttingPatchIDs_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::meshToMesh::srcMesh() const
{
    return srcMesh_;
}


inline const Foam::polyMesh& Foam::meshToMesh::tgtMesh() const
{
    return tgtMesh_;
}


inline Foam::label Foam::meshToMesh::singleProcess() const
{
    return singleProcess_;
}


inline bool Foam::meshToMesh::isSingleProcess() const
{
    return singleProcess_ != -1;
}


// ************************************************************************* //
