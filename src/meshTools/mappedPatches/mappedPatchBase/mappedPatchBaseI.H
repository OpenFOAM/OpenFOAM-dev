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

#include "mappedPatchBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline bool Foam::mappedPatchBase::nbrPatchIsMapped() const
{
    return isA<mappedPatchBase>(nbrPolyPatch());
}


inline const Foam::mappedPatchBase&
Foam::mappedPatchBase::nbrMappedPatch() const
{
    return refCast<const mappedPatchBase>(nbrPolyPatch());
}


inline bool Foam::mappedPatchBase::symmetric() const
{
    return
        nbrPatchIsMapped()
     && nbrMappedPatch().nbrRegionName_ == patch_.boundaryMesh().mesh().name()
     && nbrMappedPatch().nbrPatchName_ == patch_.name()
     && !usingTree_
     && !nbrMappedPatch().usingTree_
     && nbrMappedPatch().patchToPatchPtr_->type() == patchToPatchPtr_->type();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::mappedPatchBase::nbrRegionName() const
{
    if (coupleGroup_.valid() && nbrRegionName_.empty())
    {
        // Try and use the couple group to find nbrPatch and nbrRegion. Note:
        // This call sets nbrRegionName_.
        const label nbrPatchID =
            coupleGroup_.findOtherPatchID(patch_, nbrRegionName_);

        nbrPatchName_ = nbrMesh().boundaryMesh()[nbrPatchID].name();
    }

    return nbrRegionName_;
}


inline const Foam::word& Foam::mappedPatchBase::nbrPatchName() const
{
    if (nbrPatchName_.empty())
    {
        nbrRegionName();
    }

    return nbrPatchName_;
}


inline bool Foam::mappedPatchBase::sameRegion() const
{
    return nbrRegionName() == patch_.boundaryMesh().mesh().name();
}


inline bool Foam::mappedPatchBase::samePatch() const
{
    return sameRegion() && nbrPatchName() == patch_.name();
}


inline bool Foam::mappedPatchBase::sameUntransformedPatch() const
{
    return
        samePatch()
     && (
            !transform_.transformComplete()
         || !transform_.transform().transformsPosition()
        );
}


// ************************************************************************* //
