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

#include "volume.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneGenType>
Foam::zoneSet Foam::zoneGenerators::volume::generate
(
    const ZoneGenType& zoneGen
) const
{
    labelList indices;

    switch (zoneType_)
    {
        case zoneTypes::point:
        {
            return zoneSet
            (
                new pointZone
                (
                    zoneName_,
                    select(zoneGen, mesh_.points()),
                    mesh_.pointZones(),
                    moveUpdate_,
                    true
                )
            );
        }

        case zoneTypes::cell:
        {
            return zoneSet
            (
                new cellZone
                (
                    zoneName_,
                    select(zoneGen, mesh_.cellCentres()),
                    mesh_.cellZones(),
                    moveUpdate_,
                    true
                )
            );
        }

        case zoneTypes::face:
        {
            const labelList faceIndices(select(zoneGen, mesh_.faceCentres()));

            return zoneSet
            (
                new faceZone
                (
                    zoneName_,
                    faceIndices,
                    mesh_.faceZones(),
                    moveUpdate_,
                    true
                )
            );
        }
    }

    return zoneSet();
}


// ************************************************************************* //
