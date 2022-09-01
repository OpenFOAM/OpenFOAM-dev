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

inline bool Foam::mappedPatchBase::sampleIsMappedPatch() const
{
    switch (mode_)
    {
        case NEARESTCELL:
        case NEARESTFACE:
            return false;
        case NEARESTPATCHFACE:
        case NEARESTPATCHFACEAMI:
        case PATCHTOPATCH:
            return isA<mappedPatchBase>(samplePolyPatch());
    }
    return false;
}


inline const Foam::mappedPatchBase&
Foam::mappedPatchBase::sampleMappedPatch() const
{
    return refCast<const mappedPatchBase>(samplePolyPatch());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::mappedPatchBase::sampleMode&
Foam::mappedPatchBase::mode() const
{
    return mode_;
}


inline const Foam::word& Foam::mappedPatchBase::sampleRegion() const
{
    if (sampleRegion_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a regionName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find samplePatch and sampleRegion
        const label samplePatchID =
            coupleGroup_.findOtherPatchID(patch_, sampleRegion_);

        sameRegion_ = sampleRegion_ == patch_.boundaryMesh().mesh().name();

        samplePatch_ = sampleMesh().boundaryMesh()[samplePatchID].name();
    }

    return sampleRegion_;
}


inline const Foam::word& Foam::mappedPatchBase::samplePatch() const
{
    if (samplePatch_.empty())
    {
        sampleRegion();
    }

    return samplePatch_;
}


inline bool Foam::mappedPatchBase::sameRegion() const
{
    return sameRegion_;
}


inline const Foam::distributionMap& Foam::mappedPatchBase::map() const
{
    if (mapPtr_.empty())
    {
        calcMapping();
    }

    return mapPtr_();
}


inline const Foam::labelList& Foam::mappedPatchBase::mapIndices() const
{
    if (mapPtr_.empty())
    {
        calcMapping();
    }

    return mapIndices_;
}


// ************************************************************************* //
