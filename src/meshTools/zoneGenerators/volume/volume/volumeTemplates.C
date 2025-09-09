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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ZoneType, class UnaryOp, class ZoneGenType>
inline Foam::labelList Foam::zoneGenerators::volume::select
(
    const ZoneGenType& zoneGen,
    const vectorField& pts,
    const UnaryOp& uop
) const
{
    labelList indices(pts.size());

    label nInZone = 0;
    forAll(pts, i)
    {
        if (uop(zoneGen.contains(pts[i])))
        {
            indices[nInZone++] = i;
        }
    }

    indices.setSize(nInZone);

    return indices;
}


template<class ZoneType, class UnaryOp, class ZoneGenType>
inline Foam::labelList Foam::zoneGenerators::volume::select
(
    const ZoneGenType& zoneGen,
    const zoneGeneratorList& zoneGenerators,
    const vectorField& pts,
    const UnaryOp& uop
) const
{
    if (!zoneGenerators.size())
    {
        return select<ZoneType, UnaryOp>(zoneGen, pts, uop);
    }

    boolList selected(pts.size(), false);

    forAll(zoneGenerators, zgi)
    {
        zoneSet zs(zoneGenerators[zgi].generate());

        const labelList& zsIndices(zs.zone<ZoneType>());

        forAll(zsIndices, i)
        {
            if (uop(zoneGen.contains(pts[zsIndices[i]])))
            {
                selected[zsIndices[i]] = true;
            }
        }
    }

    return indices(selected);
}


template<class ZoneType, class ZoneGenType>
inline Foam::labelList Foam::zoneGenerators::volume::selectOp
(
    const ZoneGenType& zoneGen,
    const zoneGeneratorList& zoneGenerators,
    const vectorField& pts
) const
{
    if (select_ == selection::inside)
    {
        return zoneGen.template select<ZoneType>
        (
            zoneGen,
            zoneGenerators,
            pts,
            nopOp<bool>()
        );
    }
    else
    {
        return zoneGen.template select<ZoneType>
        (
            zoneGen,
            zoneGenerators,
            pts,
            notOp<bool>()
        );
    }
}


template<class UnaryOp, class ZoneGenType>
inline Foam::labelList Foam::zoneGenerators::volume::select
(
    const ZoneGenType& zoneGen,
    const zoneGeneratorList& zoneGenerators,
    const vectorField& pts,
    boolList& flipMap,
    const UnaryOp& uop
) const
{
    if (!zoneGenerators.size())
    {
        return select<faceZone>(zoneGen, pts, uop);
    }

    bool oriented = true;
    boolList selected(pts.size(), false);

    forAll(zoneGenerators, zgi)
    {
        zoneSet zs(zoneGenerators[zgi].generate());

        const faceZone& fZone = zs.fZone();
        const labelList& zsIndices(fZone);

        if (oriented && fZone.oriented())
        {
            flipMap.setSize(mesh_.nFaces(), false);
            const boolList& zsFlipMap(fZone.flipMap());

            forAll(zsIndices, i)
            {
                if (uop(zoneGen.contains(pts[zsIndices[i]])))
                {
                    selected[zsIndices[i]] = true;
                    flipMap[zsIndices[i]] = zsFlipMap[i];
                }
            }
        }
        else
        {
            oriented = false;
            flipMap.clear();

            forAll(zsIndices, i)
            {
                if (uop(zoneGen.contains(pts[zsIndices[i]])))
                {
                    selected[zsIndices[i]] = true;
                }
            }
        }
    }

    labelList faceIndices(indices(selected));

    if (oriented)
    {
        flipMap = boolList(flipMap, faceIndices);
    }

    return faceIndices;
}


template<class ZoneGenType>
inline Foam::labelList Foam::zoneGenerators::volume::selectOp
(
    const ZoneGenType& zoneGen,
    const zoneGeneratorList& zoneGenerators,
    const vectorField& pts,
    boolList& flipMap
) const
{
    if (select_ == selection::inside)
    {
        return zoneGen.select
        (
            zoneGen,
            zoneGenerators,
            pts,
            flipMap,
            nopOp<bool>()
        );
    }
    else
    {
        return zoneGen.select
        (
            zoneGen,
            zoneGenerators,
            pts,
            flipMap,
            notOp<bool>()
        );
    }
}


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
                    zoneGen.template selectOp<pointZone>
                    (
                        zoneGen,
                        zoneGenerators_,
                        mesh_.points()
                    ),
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
                    zoneGen.template selectOp<cellZone>
                    (
                        zoneGen,
                        zoneGenerators_,
                        mesh_.cellCentres()
                    ),
                    mesh_.cellZones(),
                    moveUpdate_,
                    true
                )
            );
        }

        case zoneTypes::face:
        {
            boolList flipMap;
            const labelList faceIndices
            (
                zoneGen.selectOp
                (
                    zoneGen,
                    zoneGenerators_,
                    mesh_.faceCentres(),
                    flipMap
                )
            );

            return zoneSet
            (
                flipMap.size()
                  ? new faceZone
                    (
                        zoneName_,
                        faceIndices,
                        flipMap,
                        mesh_.faceZones(),
                        moveUpdate_,
                        true
                    )
                  : new faceZone
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
