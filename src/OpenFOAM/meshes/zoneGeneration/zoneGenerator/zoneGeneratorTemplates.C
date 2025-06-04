/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "zoneGenerator.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ZoneGenType>
inline Foam::labelList Foam::zoneGenerator::zoneGenerator::select
(
    const ZoneGenType& zoneGen,
    const vectorField& pts
) const
{
    labelList indices(pts.size());

    label nInZone = 0;
    forAll(pts, i)
    {
        if (zoneGen.contains(pts[i]))
        {
            indices[nInZone++] = i;
        }
    }

    indices.setSize(nInZone);

    return indices;
}


// ************************************************************************* //
