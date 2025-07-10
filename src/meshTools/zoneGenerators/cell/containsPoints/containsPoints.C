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

#include "containsPoints.H"
#include "meshSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(containsPoints, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            containsPoints,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::containsPoints::containsPoints
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    points_(dict.lookup<List<point>>("points", dimLength))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::containsPoints::~containsPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::containsPoints::generate() const
{
    const meshSearch& searchEngine = meshSearch::New(mesh_);

    labelList cellIndices(points_.size());

    label czi = 0;

    forAll(points_, i)
    {
        const label celli = searchEngine.findCell(points_[i]);
        if (celli >= 0)
        {
            cellIndices[czi++] = celli;
        }

        const label globalCelli = returnReduce(celli, maxOp<label>());
        if (globalCelli < 0)
        {
            WarningInFunction
                << "Unable to find cell that contains point " << points_[i]
                << endl;
        }
    }

    cellIndices.setSize(czi);

    return zoneSet
    (
        new cellZone
        (
            zoneName_,
            cellIndices,
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
