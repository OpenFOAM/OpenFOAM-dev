/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "sampledIsoSurfaceSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(sampledIsoSurfaceSurface, 0);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::cutPolyIsoSurface>
Foam::sampledSurfaces::sampledIsoSurfaceSurface::calcIsoSurf
(
    const scalarField& pointValues,
    const scalar isoValue
) const
{
    if (zone_.all())
    {
        return autoPtr<cutPolyIsoSurface>
        (
            new cutPolyIsoSurface(mesh(), pointValues, isoValue)
        );
    }
    else
    {
        return autoPtr<cutPolyIsoSurface>
        (
            new cutPolyIsoSurface(mesh(), pointValues, isoValue, zone_.zone())
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledIsoSurfaceSurface::sampledIsoSurfaceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    zone_(mesh, dict, true),
    isoSurfPtr_(nullptr),
    isoSurfTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledIsoSurfaceSurface::~sampledIsoSurfaceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::sampledIsoSurfaceSurface::update() const
{
    // Quick return if no update needed
    if (!needsUpdate())
    {
        return false;
    }

    // Clear any information in the base class
    sampledSurface::clearGeom();

    zone_.regenerate();

    // Update the iso surface
    isoSurfPtr_.reset(calcIsoSurf().ptr());

    // Set the time index
    isoSurfTimeIndex_ = mesh().time().timeIndex();

    return true;
}


bool Foam::sampledSurfaces::sampledIsoSurfaceSurface::update()
{
    return static_cast<const sampledIsoSurfaceSurface&>(*this).update();
}


#define IMPLEMENT_SAMPLE(Type, nullArg)                                        \
    Foam::tmp<Foam::Field<Foam::Type>>                                         \
    Foam::sampledSurfaces::sampledIsoSurfaceSurface::sample                    \
    (                                                                          \
        const VolField<Type>& vField                                           \
    ) const                                                                    \
    {                                                                          \
        return sampleField(vField);                                            \
    }
FOR_ALL_FIELD_TYPES(IMPLEMENT_SAMPLE);
#undef IMPLEMENT_SAMPLE


#define IMPLEMENT_INTERPOLATE(Type, nullArg)                                   \
    Foam::tmp<Foam::Field<Foam::Type>>                                         \
    Foam::sampledSurfaces::sampledIsoSurfaceSurface::interpolate               \
    (                                                                          \
        const interpolation<Type>& interpolator                                \
    ) const                                                                    \
    {                                                                          \
        return interpolateField(interpolator);                                 \
    }
FOR_ALL_FIELD_TYPES(IMPLEMENT_INTERPOLATE);
#undef IMPLEMENT_INTERPOLATE


void Foam::sampledSurfaces::sampledIsoSurfaceSurface::movePoints()
{
    sampledSurface::clearGeom();
    zone_.movePoints();
    isoSurfPtr_.clear();
    isoSurfTimeIndex_ = -1;
}


void Foam::sampledSurfaces::sampledIsoSurfaceSurface::topoChange
(
    const polyTopoChangeMap& map
)
{
    sampledSurface::clearGeom();
    zone_.topoChange(map);
    isoSurfPtr_.clear();
    isoSurfTimeIndex_ = -1;
}


void Foam::sampledSurfaces::sampledIsoSurfaceSurface::mapMesh
(
    const polyMeshMap& map
)
{
    sampledSurface::clearGeom();
    zone_.mapMesh(map);
    isoSurfPtr_.clear();
    isoSurfTimeIndex_ = -1;
}


void Foam::sampledSurfaces::sampledIsoSurfaceSurface::distribute
(
    const polyDistributionMap& map
)
{
    sampledSurface::clearGeom();
    zone_.distribute(map);
    isoSurfPtr_.clear();
    isoSurfTimeIndex_ = -1;
}


// ************************************************************************* //
