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

#include "orientZoneGenerator.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(orientZoneGenerator, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            orientZoneGenerator,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::orientZoneGenerator::orientZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerator_(zoneGenerator::New(mesh, dict)),
    normal_(dict.lookup<vector>("normal", dimless))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::orientZoneGenerator::~orientZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::orientZoneGenerator::generate() const
{
    const labelList faceIndices(zoneGenerator_->generate().fZone());
    boolList flipMap(faceIndices.size(), false);

    const vectorField& faceAreas = mesh_.faceAreas();

    forAll(faceIndices, fi)
    {
        const label facei = faceIndices[fi];

        if ((faceAreas[facei] & normal_) < 0)
        {
            flipMap[fi] = true;
        }
    }

    return zoneSet
    (
        new faceZone
        (
            zoneName_,
            faceIndices,
            flipMap,
            mesh_.faceZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
