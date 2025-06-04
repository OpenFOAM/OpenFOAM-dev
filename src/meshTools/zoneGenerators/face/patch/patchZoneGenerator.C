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

#include "patchZoneGenerator.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(patchZoneGenerator, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            patchZoneGenerator,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::patchZoneGenerator::patchZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    patchSet_(mesh.boundaryMesh().patchSet(dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::patchZoneGenerator::~patchZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::patchZoneGenerator::generate() const
{
    boolList selectedFaces(mesh_.nFaces(), false);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label patchi = iter.key();
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        for
        (
            label facei = pp.start();
            facei < pp.start() + pp.size();
            facei++
        )
        {
            selectedFaces[facei] = true;
        }
    }

    const labelList faceIndices(indices(selectedFaces));

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            faceIndices,
            boolList(faceIndices.size(), false),
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
