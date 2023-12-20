/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

inline const Foam::cellsToCells& Foam::meshToMesh::cellsInterpolation() const
{
    return cellsInterpolation_();
}


inline const Foam::cellsToCellsStabilisation&
Foam::meshToMesh::srcCellsStabilisation() const
{
    if (!srcCellsStabilisation_.valid())
    {
        srcCellsStabilisation_.reset(new cellsToCellsStabilisation());
        srcCellsStabilisation_->update
        (
            srcMesh_,
            cellsInterpolation_->srcCoupled()
        );
    }

    return srcCellsStabilisation_();
}


inline const Foam::cellsToCellsStabilisation&
Foam::meshToMesh::tgtCellsStabilisation() const
{
    if (!tgtCellsStabilisation_.valid())
    {
        tgtCellsStabilisation_.reset(new cellsToCellsStabilisation());
        tgtCellsStabilisation_->update
        (
            tgtMesh_,
            cellsInterpolation_->tgtCoupled()
        );
    }

    return tgtCellsStabilisation_();
}


inline const Foam::List<Foam::labelPair>& Foam::meshToMesh::patchIndices() const
{
    return patchIndices_;
}


inline const Foam::patchToPatch& Foam::meshToMesh::patchInterpolation
(
    const label i
) const
{
    return patchInterpolations_[i];
}


inline const Foam::patchToPatchStabilisation&
Foam::meshToMesh::srcPatchStabilisation(const label i) const
{
    if (!srcPatchStabilisations_.set(i))
    {
        const label srcPatchi = patchIndices_[i].first();
        const polyPatch& srcPp = srcMesh_.boundaryMesh()[srcPatchi];

        srcPatchStabilisations_.set(i, new patchToPatchStabilisation());
        srcPatchStabilisations_[i].update
        (
            srcPp,
            patchInterpolations_[i].srcCoupled()
        );
    }

    return srcPatchStabilisations_[i];
}


inline const Foam::patchToPatchStabilisation&
Foam::meshToMesh::tgtPatchStabilisation(const label i) const
{
    if (!tgtPatchStabilisations_.set(i))
    {
        const label tgtPatchi = patchIndices_[i].first();
        const polyPatch& tgtPp = tgtMesh_.boundaryMesh()[tgtPatchi];

        tgtPatchStabilisations_.set(i, new patchToPatchStabilisation());
        tgtPatchStabilisations_[i].update
        (
            tgtPp,
            patchInterpolations_[i].tgtCoupled()
        );
    }

    return tgtPatchStabilisations_[i];
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


// ************************************************************************* //
