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

#include "mappedPatchBaseBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::mappedPatchBaseBase::nbrRegionName() const
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


inline const Foam::word& Foam::mappedPatchBaseBase::nbrPatchName() const
{
    if (nbrPatchName_.empty())
    {
        nbrRegionName();
    }

    return nbrPatchName_;
}


inline const Foam::transformer& Foam::mappedPatchBaseBase::transform() const
{
    return transform_.transform();
}


inline bool Foam::mappedPatchBaseBase::sameRegion() const
{
    return nbrRegionName() == patch_.boundaryMesh().mesh().name();
}


inline bool Foam::mappedPatchBaseBase::samePatch() const
{
    return sameRegion() && nbrPatchName() == patch_.name();
}


inline bool Foam::mappedPatchBaseBase::sameUntransformedPatch() const
{
    return
        samePatch()
     && (
            !transform_.transformComplete()
         || !transform_.transform().transformsPosition()
        );
}


// ************************************************************************* //
