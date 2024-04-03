/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "faceZoneList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneList, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::boolList Foam::faceZoneList::zonesFlipFace
(
    const label facei,
    const labelList& faceiZones
)
{
    labelList zones(whichZones(facei));
    boolList flipFaces(zones.size());

    forAll(zones, zi)
    {
        const faceZone& fz = this->operator[](zi);
        flipFaces[zi] = fz.flipMap()[fz.localIndex(facei)];
    }

    return flipFaces;
}


void Foam::faceZoneList::insert(const List<Map<bool>>& zonesIndices)
{
    PtrList<faceZone>& zones = *this;

    if (zonesIndices.size() != zones.size())
    {
        FatalErrorInFunction
            << "zonesIndices.size() " << zonesIndices.size()
            << " != number of zones " << zones.size()
            << exit(FatalError);
    }

    forAll(zonesIndices, zonei)
    {
        zones[zonei].insert(zonesIndices[zonei]);
    }
}


// ************************************************************************* //
