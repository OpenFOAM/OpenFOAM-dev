/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "surfaceZonesInfo.H"
#include "searchableSurface.H"
#include "searchableSurfaces.H"
#include "polyMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::surfaceZonesInfo::areaSelectionAlgo,
        4
    >::names[] =
    {
        "inside",
        "outside",
        "insidePoint",
        "none"
    };
}
const Foam::NamedEnum<Foam::surfaceZonesInfo::areaSelectionAlgo, 4>
    Foam::surfaceZonesInfo::areaSelectionAlgoNames;


namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::surfaceZonesInfo::faceZoneType,
        3
    >::names[] =
    {
        "internal",
        "baffle",
        "boundary"
    };
}
const Foam::NamedEnum<Foam::surfaceZonesInfo::faceZoneType, 3>
    Foam::surfaceZonesInfo::faceZoneTypeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceZonesInfo::surfaceZonesInfo
(
    const searchableSurface& surface,
    const dictionary& surfacesDict
)
:
    faceZoneName_(),
    cellZoneName_(),
    zoneInside_(NONE),
    zoneInsidePoint_(point::min),
    faceType_(INTERNAL)
{
    // Global zone names per surface
    if (surfacesDict.readIfPresent("faceZone", faceZoneName_))
    {
        // Read optional entry to determine inside of faceZone

        word method;
        bool hasSide =
            surfacesDict.readIfPresent("mode", method)
         || surfacesDict.readIfPresent("cellZoneInside", method);

        if (hasSide)
        {
            zoneInside_ = areaSelectionAlgoNames[method];
            if (zoneInside_ == INSIDEPOINT)
            {
                surfacesDict.lookup("insidePoint") >> zoneInsidePoint_;
            }
        }

        // Read optional cellZone name

        if (surfacesDict.readIfPresent("cellZone", cellZoneName_))
        {
            if
            (
                (
                    zoneInside_ == INSIDE
                 || zoneInside_ == OUTSIDE
                )
            && !surface.hasVolumeType()
            )
            {
                IOWarningInFunction
                (
                    surfacesDict
                )   << "Illegal entry zoneInside "
                    << areaSelectionAlgoNames[zoneInside_]
                    << " for faceZone "
                    << faceZoneName_
                    << " since surface is not closed." << endl;
            }
        }
        else if (hasSide)
        {
            IOWarningInFunction
            (
                surfacesDict
            )   << "Unused entry zoneInside for faceZone "
                << faceZoneName_
                << " since no cellZone specified."
                << endl;
        }

        // How to handle faces on faceZone
        word faceTypeMethod;
        if (surfacesDict.readIfPresent("faceType", faceTypeMethod))
        {
            faceType_ = faceZoneTypeNames[faceTypeMethod];
        }
    }
}


Foam::surfaceZonesInfo::surfaceZonesInfo
(
    const word& faceZoneName,
    const word& cellZoneName,
    const areaSelectionAlgo& zoneInside,
    const point& zoneInsidePoint,
    const faceZoneType& faceType
)
:
    faceZoneName_(faceZoneName),
    cellZoneName_(cellZoneName),
    zoneInside_(zoneInside),
    zoneInsidePoint_(zoneInsidePoint),
    faceType_(faceType)
{}


Foam::surfaceZonesInfo::surfaceZonesInfo(const surfaceZonesInfo& surfZone)
:
    faceZoneName_(surfZone.faceZoneName()),
    cellZoneName_(surfZone.cellZoneName()),
    zoneInside_(surfZone.zoneInside()),
    zoneInsidePoint_(surfZone.zoneInsidePoint()),
    faceType_(surfZone.faceType())
{}


Foam::labelList Foam::surfaceZonesInfo::getUnnamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
    labelList anonymousSurfaces(surfList.size());

    label i = 0;
    forAll(surfList, surfI)
    {
        if (surfList[surfI].faceZoneName().empty())
        {
            anonymousSurfaces[i++] = surfI;
        }
    }
    anonymousSurfaces.setSize(i);

    return anonymousSurfaces;
}


Foam::labelList Foam::surfaceZonesInfo::getNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
   labelList namedSurfaces(surfList.size());

    label namedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].faceZoneName().size()
        )
        {
            namedSurfaces[namedI++] = surfI;
        }
    }
    namedSurfaces.setSize(namedI);

    return namedSurfaces;
}


Foam::labelList Foam::surfaceZonesInfo::getClosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && (
                surfList[surfI].zoneInside() == surfaceZonesInfo::INSIDE
             || surfList[surfI].zoneInside() == surfaceZonesInfo::OUTSIDE
            )
         && allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::labelList Foam::surfaceZonesInfo::getUnclosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList unclosed(surfList.size());

    label unclosedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && !allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            unclosed[unclosedI++] = surfI;
        }
    }
    unclosed.setSize(unclosedI);

    return unclosed;
}


Foam::labelList Foam::surfaceZonesInfo::getAllClosedNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList,
    const searchableSurfaces& allGeometry,
    const labelList& surfaces
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && allGeometry[surfaces[surfI]].hasVolumeType()
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::labelList Foam::surfaceZonesInfo::getInsidePointNamedSurfaces
(
    const PtrList<surfaceZonesInfo>& surfList
)
{
    labelList closed(surfList.size());

    label closedI = 0;
    forAll(surfList, surfI)
    {
        if
        (
            surfList.set(surfI)
         && surfList[surfI].cellZoneName().size()
         && surfList[surfI].zoneInside() == surfaceZonesInfo::INSIDEPOINT
        )
        {
            closed[closedI++] = surfI;
        }
    }
    closed.setSize(closedI);

    return closed;
}


Foam::labelList Foam::surfaceZonesInfo::addCellZonesToMesh
(
    const PtrList<surfaceZonesInfo>& surfList,
    const labelList& namedSurfaces,
    polyMesh& mesh
)
{
    labelList surfaceToCellZone(surfList.size(), -1);

    cellZoneMesh& cellZones = mesh.cellZones();

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        const word& cellZoneName = surfList[surfI].cellZoneName();

        if (cellZoneName != word::null)
        {
            label zoneI = cellZones.findZoneID(cellZoneName);

            if (zoneI == -1)
            {
                zoneI = cellZones.size();
                cellZones.setSize(zoneI+1);
                cellZones.set
                (
                    zoneI,
                    new cellZone
                    (
                        cellZoneName,   // name
                        labelList(0),   // addressing
                        zoneI,          // index
                        cellZones       // cellZoneMesh
                    )
                );
            }

            surfaceToCellZone[surfI] = zoneI;
        }
    }

    // Check they are synced
    List<wordList> allCellZones(Pstream::nProcs());
    allCellZones[Pstream::myProcNo()] = cellZones.names();
    Pstream::gatherList(allCellZones);
    Pstream::scatterList(allCellZones);

    for (label proci = 1; proci < allCellZones.size(); proci++)
    {
        if (allCellZones[proci] != allCellZones[0])
        {
            FatalErrorInFunction
                << "Zones not synchronised among processors." << nl
                << " Processor0 has cellZones:" << allCellZones[0]
                << " , processor" << proci
                << " has cellZones:" << allCellZones[proci]
                << exit(FatalError);
        }
    }

    return surfaceToCellZone;
}


Foam::labelList Foam::surfaceZonesInfo::addFaceZonesToMesh
(
    const PtrList<surfaceZonesInfo>& surfList,
    const labelList& namedSurfaces,
    polyMesh& mesh
)
{
    labelList surfaceToFaceZone(surfList.size(), -1);

    faceZoneMesh& faceZones = mesh.faceZones();

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];

        const word& faceZoneName = surfList[surfI].faceZoneName();

        label zoneI = faceZones.findZoneID(faceZoneName);

        if (zoneI == -1)
        {
            zoneI = faceZones.size();
            faceZones.setSize(zoneI+1);
            faceZones.set
            (
                zoneI,
                new faceZone
                (
                    faceZoneName,   // name
                    labelList(0),   // addressing
                    boolList(0),    // flipmap
                    zoneI,          // index
                    faceZones       // faceZoneMesh
                )
            );
        }

        surfaceToFaceZone[surfI] = zoneI;
    }

    // Check they are synced
    List<wordList> allFaceZones(Pstream::nProcs());
    allFaceZones[Pstream::myProcNo()] = faceZones.names();
    Pstream::gatherList(allFaceZones);
    Pstream::scatterList(allFaceZones);

    for (label proci = 1; proci < allFaceZones.size(); proci++)
    {
        if (allFaceZones[proci] != allFaceZones[0])
        {
            FatalErrorInFunction
                << "Zones not synchronised among processors." << nl
                << " Processor0 has faceZones:" << allFaceZones[0]
                << " , processor" << proci
                << " has faceZones:" << allFaceZones[proci]
                << exit(FatalError);
        }
    }

    return surfaceToFaceZone;
}


// ************************************************************************* //
