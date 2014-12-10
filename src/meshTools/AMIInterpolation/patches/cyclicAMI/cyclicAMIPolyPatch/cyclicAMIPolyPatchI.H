/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::word& Foam::cyclicAMIPolyPatch::neighbPatchName() const
{
    if (nbrPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        nbrPatchName_ = boundaryMesh()[patchID].name();
    }
    return nbrPatchName_;
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::rotationAxis() const
{
    return rotationAxis_;
}


inline const Foam::point& Foam::cyclicAMIPolyPatch::rotationCentre() const
{
    return rotationCentre_;
}


inline const Foam::vector& Foam::cyclicAMIPolyPatch::separationVector() const
{
    return separationVector_;
}


// ************************************************************************* //
