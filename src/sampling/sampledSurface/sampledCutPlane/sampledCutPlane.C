/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "sampledCutPlane.H"
#include "emptyFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(cutPlane, 0);
    addToRunTimeSelectionTable(sampledSurface, cutPlane, word);

    // Backwards compatible lookup as "cuttingPlane"
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        cutPlane,
        word,
        cuttingPlane
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::cutPolyIsoSurface>
Foam::sampledSurfaces::cutPlane::calcIsoSurf() const
{
    // Compute the distance from the mesh points to the plane
    scalarField pointDistance(mesh().nPoints());
    {
        forAll(mesh().points(), pointi)
        {
            pointDistance[pointi] =
                plane_.signedDistance(mesh().points()[pointi]);
        }
    }

    // Construct an iso-surface at the given distance
    return autoPtr<cutPolyIsoSurface>
    (
        new cutPolyIsoSurface(mesh(), pointDistance, 0, zoneIDs())
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::cutPlane::cutPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledIsoSurfaceSurface(name, mesh, dict),
    plane_(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::cutPlane::~cutPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::cutPlane::needsUpdate() const
{
    return timeIndex() == -1;
}


void Foam::sampledSurfaces::cutPlane::print(Ostream& os) const
{
    os  << "cutPlane: " << name() << " :"
        << "  plane:" << plane_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
