/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "ListOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::labelList Foam::emptyLabelList = Foam::labelList(0);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::invert
(
    const label len,
    const labelUList& map
)
{
    labelList inverse(len, -1);

    forAll(map, i)
    {
        label newPos = map[i];

        if (newPos >= 0)
        {
            if (inverse[newPos] >= 0)
            {
                FatalErrorInFunction
                    << "Map is not one-to-one. At index " << i
                    << " element " << newPos << " has already occurred before"
                    << nl << "Please use invertOneToMany instead"
                    << abort(FatalError);
            }

            inverse[newPos] = i;
        }
    }
    return inverse;
}


Foam::labelListList Foam::invertOneToMany
(
    const label len,
    const labelUList& map
)
{
    labelList nElems(len, 0);

    forAll(map, i)
    {
        if (map[i] >= 0)
        {
            nElems[map[i]]++;
        }
    }

    labelListList inverse(len);

    forAll(nElems, i)
    {
        inverse[i].setSize(nElems[i]);
        nElems[i] = 0;
    }

    forAll(map, i)
    {
        label newI = map[i];

        if (newI >= 0)
        {
            inverse[newI][nElems[newI]++] = i;
        }
    }

    return inverse;
}


Foam::labelList Foam::identity(const label len)
{
    labelList map(len);

    forAll(map, i)
    {
        map[i] = i;
    }
    return map;
}


// ************************************************************************* //
