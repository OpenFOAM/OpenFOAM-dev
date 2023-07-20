/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "meshingSurface.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<meshingSurface::surfaceType, 5>::names[] =
        {"wall", "external", "cellZone", "rotatingZone", "baffle"};
}


const Foam::NamedEnum<Foam::meshingSurface::surfaceType, 5>
    Foam::meshingSurface::surfaceTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshingSurface::nSurfaceParts(const triSurfaceMesh& surf)
{
    labelList faceZone;
    return surf.markZones(boolList(surf.nEdges(), false), faceZone);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshingSurface::meshingSurface()
:
    path_(),
    file_(),
    name_(),
    type_(surfaceType::wall),
    boundBox_(),
    closed_(),
    inletRegions_(),
    outletRegions_()
{}


Foam::meshingSurface::meshingSurface(const fileName& file, const Time& time)
:
    path_(file.path()),
    file_(file.name()),
    name_(word(file_.lessExt())),
    type_(surfaceType::wall),
    boundBox_(),
    closed_(),
    nParts_(),
    regions_(),
    inletRegions_(),
    outletRegions_()
{
    triSurfaceMesh surf
    (
        IOobject
        (
            file_,
            path_,
            time
        ),
        triSurface(path_/file_)
    );

    surf.cleanup(false);

    boundBox_ = Foam::boundBox(surf.points(), false);
    closed_ = surf.hasVolumeType();
    nParts_ = nSurfaceParts(surf);

    forAll(surf.patches(), i)
    {
        regions_.append(surf.patches()[i].name());
    }

    // Add any region with name beginning "inlet" to inletRegions_
    forAll(regions_, r)
    {
        if (!strncmp(regions_[r].c_str(), "inlet", 5))
        {
            inletRegions_.append(regions_[r]);
        }
    }

    // Add any region with name beginning "outlet" to outletRegions_
    forAll(regions_, r)
    {
        if (!strncmp(regions_[r].c_str(), "outlet", 6))
        {
            outletRegions_.append(regions_[r]);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshingSurface::~meshingSurface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::meshingSurface::isSurfaceExt(const Foam::fileName& file)
{
    const word ext(file.ext());

    if (ext == "stl" || ext == "stlb" || ext == "obj" || ext == "vtk")
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
