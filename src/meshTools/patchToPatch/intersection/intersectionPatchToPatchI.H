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

#include "intersectionPatchToPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const
Foam::List<Foam::DynamicList<Foam::patchToPatches::intersection::couple>>&
Foam::patchToPatches::intersection::srcCouples() const
{
    return srcCouples_;
}


inline const Foam::List<Foam::patchToPatches::intersection::part>&
Foam::patchToPatches::intersection::srcEdgeParts() const
{
    return srcEdgeParts_;
}


inline const Foam::List<Foam::patchToPatches::intersection::part>&
Foam::patchToPatches::intersection::srcErrorParts() const
{
    return srcErrorParts_;
}


inline const
Foam::List<Foam::DynamicList<Foam::patchToPatches::intersection::couple>>&
Foam::patchToPatches::intersection::tgtCouples() const
{
    return tgtCouples_;
}


inline const Foam::List<Foam::scalar>&
Foam::patchToPatches::intersection::srcCoverage() const
{
    return srcCoverage_;
}


inline const Foam::List<Foam::scalar>&
Foam::patchToPatches::intersection::tgtCoverage() const
{
    return tgtCoverage_;
}


// ************************************************************************* //
