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

Foam::labelList Foam::zoneGenerators::insideSurface::select
(
    const insideSurface& zoneGen,
    const vectorField& pts
) const
{
    labelList indices(pts.size());

    List<volumeType> volType;
    surfacePtr_->getVolumeType(pts, volType);

    label nInZone = 0;
    forAll(volType, i)
    {
        if (volType[i] == volumeType::inside)
        {
            indices[nInZone++] = i;
        }
    }

    indices.setSize(nInZone);

    return indices;
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
