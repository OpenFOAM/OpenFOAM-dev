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

#include "insideSurface.H"
#include "Time.H"
#include "polyMesh.H"
#include "searchableSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(insideSurface, 0);
        addToRunTimeSelectionTable(zoneGenerator, insideSurface, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class UnaryOp>
Foam::labelList Foam::zoneGenerators::insideSurface::select
(
    const insideSurface& zoneGen,
    const vectorField& pts,
    const UnaryOp& uop
) const
{
    labelList indices(pts.size());

    List<volumeType> volType;
    surfacePtr_->getVolumeType(pts, volType);

    label nInZone = 0;
    forAll(volType, i)
    {
        if (uop(volType[i] == volumeType::inside))
        {
            indices[nInZone++] = i;
        }
    }

    indices.setSize(nInZone);

    return indices;
}


template<class ZoneType, class UnaryOp>
Foam::labelList Foam::zoneGenerators::insideSurface::select
(
    const insideSurface& zoneGen,
    const zoneGeneratorList& zoneGenerators,
    const vectorField& pts,
    const UnaryOp& uop
) const
{
    if (!zoneGenerators.size())
    {
        return zoneGen.select<ZoneType>(zoneGen, pts, uop);
    }

    boolList selected(pts.size(), false);

    forAll(zoneGenerators, zgi)
    {
        zoneSet zs(zoneGenerators[zgi].generate());

        const labelList& zsIndices(zs.zone<ZoneType>());

        const pointField zonePoints(pts, zsIndices);
        List<volumeType> volType;
        surfacePtr_->getVolumeType(zonePoints, volType);

        forAll(zsIndices, i)
        {
            if (uop(volType[i] == volumeType::inside))
            {
                selected[zsIndices[i]] = true;
            }
        }
    }

    return indices(selected);
}


template<class UnaryOp>
Foam::labelList Foam::zoneGenerators::insideSurface::select
(
    const insideSurface& zoneGen,
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

        const pointField zonePoints(pts, zsIndices);
        List<volumeType> volType;
        surfacePtr_->getVolumeType(zonePoints, volType);

        if (oriented && fZone.oriented())
        {
            flipMap.setSize(mesh_.nFaces(), false);
            const boolList& zsFlipMap(fZone.flipMap());

            forAll(zsIndices, i)
            {
                if (uop(volType[i] == volumeType::inside))
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
                if (uop(volType[i] == volumeType::inside))
                {
                    selected[zsIndices[i]] = true;
                }
            }
        }
    }

    return indices(selected);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::insideSurface::insideSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    volume(name, mesh, dict),
    surfacePtr_
    (
        searchableSurface::New
        (
            word(dict.lookup("surface")),
            IOobject
            (
                dict.lookupOrDefault
                (
                    "surfaceName",
                    mesh.objectRegistry::db().name()
                ),
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::insideSurface::~insideSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::insideSurface::generate() const
{
    return volume::generate(*this);
}


// ************************************************************************* //
