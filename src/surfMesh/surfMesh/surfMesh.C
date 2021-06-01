/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "surfMesh.H"
#include "MeshedSurfaceProxy.H"

#include "Time.H"
#include "OSspecific.H"
#include "MeshedSurface.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMesh, 0);
}


Foam::word Foam::surfMesh::meshSubDir = "surfMesh";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void Foam::surfMesh::oneZone()
// {
//     word zoneName;
//
//     surfZoneList& zones = Allocator::storedIOZones();
//     if (zones.size())
//     {
//         zoneName = zones[0].name();
//     }
//     if (zoneName.empty())
//     {
//         zoneName = "zone0";
//     }
//
//     // set single default zone
//     zones.setSize(1);
//     zones[0] = surfZone
//     (
//         zoneName,
//         nFaces(),       // zone size
//         0,              // zone start
//         0               // zone index
//     );
// }


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::updatePointsRef()
{
    // Assign the reference to the points (this is truly ugly)
    reinterpret_cast<SubField<point>&>
    (
        const_cast<Field<point>&>(MeshReference::points())
    ) = reinterpret_cast<SubField<point>&>(this->storedPoints());
}


void Foam::surfMesh::updateFacesRef()
{
    // Assign the reference to the faces
    shallowCopy(this->storedFaces());
}


void Foam::surfMesh::updateRefs()
{
    this->updatePointsRef();
    this->updateFacesRef();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMesh::surfMesh(const IOobject& io, const word& surfName)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        IOobject
        (
            "surfZones",
            time().findInstance(meshDir(), "surfZones"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints())
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    pointField&& pointLst,
    faceList&& faceLst,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        move(pointLst),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        move(faceLst),
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        surfZoneList()
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints())
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    MeshedSurface<face>&& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointField(),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faceList(),
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        surfZoneList()
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints())
{
    if (debug)
    {
        Info<<"IOobject: " << io.path() << nl
            <<" name: " << io.name()
            <<" instance: " << io.instance()
            <<" local: " << io.local()
            <<" dbDir: " << io.db().dbDir() << endl;
        Info<<"creating surfMesh at instance " << instance() << endl;
        Info<<"timeName: " << instance() << endl;
    }

    // We can also send null just to initialise without allocating
    if (notNull(surf))
    {
        transfer(surf);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMesh::~surfMesh()
{
    // clearOut();
    // resetMotion();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::resetPrimitives
(
    pointField&& points,
    faceList&& faces,
    surfZoneList&& zones,
    const bool validate
)
{
    // Clear addressing.
    MeshReference::clearGeom();

    Allocator::reset(move(points), move(faces), move(zones));
    this->updateRefs();

    if (validate)
    {
        checkZones();
    }
}


void Foam::surfMesh::transfer
(
    MeshedSurface<face>& surf
)
{
    // Clear addressing.
    MeshReference::clearGeom();

    this->storedIOPoints().transfer(surf.storedPoints());
    this->storedIOFaces().transfer(surf.storedFaces());
    this->storedIOZones().transfer(surf.storedZones());

    this->updateRefs();
}


Foam::fileName Foam::surfMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::surfMesh::pointsInstance() const
{
    return this->storedIOPoints().instance();
}


const Foam::fileName& Foam::surfMesh::facesInstance() const
{
    return this->storedIOFaces().instance();
}


Foam::label Foam::surfMesh::nPoints() const
{
    return this->points().size();
}


Foam::label Foam::surfMesh::nFaces() const
{
    return this->faces().size();
}


const Foam::pointField& Foam::surfMesh::points() const
{
    return this->storedIOPoints();
}


const Foam::faceList& Foam::surfMesh::faces() const
{
    return this->storedIOFaces();
}


void Foam::surfMesh::checkZones()
{
    // extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently
    surfZoneList& zones = Allocator::storedIOZones();

    if (zones.size() <= 1)
    {
        removeZones();
    }
    else
    {
        label count = 0;
        forAll(zones, zoneI)
        {
            zones[zoneI].start() = count;
            count += zones[zoneI].size();
        }

        if (count < nFaces())
        {
            WarningInFunction
                << "more faces " << nFaces() << " than zones " << count
                << " ... extending final zone"
                << endl;

            zones.last().size() += count - nFaces();
        }
        else if (count > size())
        {
            FatalErrorInFunction
                << "more zones " << count << " than faces " << nFaces()
                << exit(FatalError);
        }
    }
}


// Add boundary patches. Constructor helper
void Foam::surfMesh::addZones
(
    const surfZoneList& srfZones,
    const bool validate
)
{
    surfZoneList& zones = Allocator::storedIOZones();

    forAll(zones, zoneI)
    {
        zones[zoneI] = surfZone(srfZones[zoneI], zoneI);
    }

    if (validate)
    {
        checkZones();
    }
}


// Remove all files and some subdirs (eg, sets)
void Foam::surfMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = db().path()/instanceDir/meshSubDir;

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"surfZones");
}

void Foam::surfMesh::removeFiles() const
{
    removeFiles(instance());
}


void Foam::surfMesh::write(const fileName& name, const surfMesh& surf)
{
    MeshedSurfaceProxy<face>
    (
        surf.points(),
        surf.faces(),
        surf.surfZones()
    ).write(name);
}


void Foam::surfMesh::write(const fileName& name)
{
    write(name, *this);
}


// ************************************************************************* //
