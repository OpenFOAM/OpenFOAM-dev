/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
     && patchToPatchIsUsed_
     && nbrMappedPatch().patchToPatchIsUsed_
     && nbrMappedPatch().patchToPatchPtr_->type() == patchToPatchPtr_->type();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::mappedPatchBase::nbrRegionName() const
{
    if (nbrRegionName_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a regionName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find nbrPatch and nbrRegion
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


// ************************************************************************* //
