/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "sampledDistanceSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::cutPolyIsoSurface>
Foam::sampledSurfaces::distanceSurface::calcIsoSurf() const
{
    // Compute the distance from the mesh points to the surface
    scalarField pointDistance(mesh().nPoints());
    {
        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            mesh().points(),
            scalarField(mesh().points().size(), great),
            nearest
        );

        if (signed_)
        {
            List<volumeType> volType;
            surfPtr_().getVolumeType(mesh().points(), volType);

            forAll(nearest, i)
            {
                if
                (
                    volType[i] != volumeType::outside
                 && volType[i] != volumeType::inside
                )
                {
                    FatalErrorInFunction
                        << "Point " << mesh().points()[i] << " could not be "
                        << "classified as either inside or outside the surface "
                        << surfPtr_->name() << exit(FatalError);
                }

                pointDistance[i] =
                    (volType[i] == volumeType::outside ? +1 : -1)
                   *mag(mesh().points()[i] - nearest[i].hitPoint());
            }
        }
        else
        {
            forAll(nearest, i)
            {
                pointDistance[i] =
                    mag(mesh().points()[i] - nearest[i].hitPoint());
            }
        }
    }

    // Construct an iso-surface at the given distance
    return autoPtr<cutPolyIsoSurface>
    (
        new cutPolyIsoSurface(mesh(), pointDistance, distance_, zoneName())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::distanceSurface::distanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledIsoSurfaceSurface(name, mesh, dict),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.lookup("surfaceType"),
            IOobject
            (
                dict.lookupOrDefault("surfaceName", name),
                mesh.time().constant(),
                searchableSurface::geometryDir(mesh.time()),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(dict.lookup<scalar>("distance")),
    signed_(readBool(dict.lookup("signed")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::distanceSurface::~distanceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::distanceSurface::needsUpdate() const
{
    return timeIndex() == -1;
}


void Foam::sampledSurfaces::distanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfPtr_().name()
        << "  distance:" << distance_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
