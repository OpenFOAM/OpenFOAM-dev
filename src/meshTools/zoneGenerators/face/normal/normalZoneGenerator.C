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

#include "normalZoneGenerator.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(normalZoneGenerator, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            normalZoneGenerator,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::normalZoneGenerator::normalZoneGenerator
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerator_(zoneGenerator::New(mesh, dict)),
    normal_(normalised(dict.lookup<vector>("normal", dimless))),
    tol_(dict.lookup<scalar>("tol", dimless))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::normalZoneGenerator::~normalZoneGenerator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::normalZoneGenerator::generate() const
{
    labelList faceIndices(zoneGenerator_->generate().fZone());
    const vectorField& faceAreas = mesh_.faceAreas();

    label fj = 0;
    forAll(faceIndices, fi)
    {
        const label facei = faceIndices[fi];

        const vector n(normalised(faceAreas[facei]));

        if (mag(1 - (n & normal_)) < tol_)
        {
            faceIndices[fj++] = facei;
        }
    }

    faceIndices.setSize(fj);

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
