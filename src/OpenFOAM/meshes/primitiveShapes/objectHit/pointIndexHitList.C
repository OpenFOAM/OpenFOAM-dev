/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "pointIndexHitList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::minDist(const List<pointIndexHit>& hitList)
{
    scalar minDist = GREAT;

    for (label hi1=0; hi1<hitList.size() - 1; ++hi1)
    {
        const pointIndexHit& pHit1 = hitList[hi1];

        if (pHit1.hit())
        {
            for (label hi2=hi1 + 1; hi2<hitList.size(); ++hi2)
            {
                const pointIndexHit& pHit2 = hitList[hi2];

                if (pHit2.hit())
                {
                    minDist =
                        min(mag(pHit1.hitPoint() - pHit2.hitPoint()), minDist);
                }
            }
        }
    }

    return minDist;
}


// ************************************************************************* //
