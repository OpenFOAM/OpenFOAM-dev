/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "cellShapes_zoneGenerator.H"
#include "polyMesh.H"
#include "hexMatcher.H"
#include "cellFeatures.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(cellShapes, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            cellShapes,
            dictionary
        );
    }
}

Foam::scalar Foam::zoneGenerators::cellShapes::featureCos
(
    Foam::cos(degToRad(10.0))
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::cellShapes::cellShapes
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    shapeTypes_(readList<word>(dict.lookup("shapes"))),
    cellModels_(shapeTypes_.size(), nullptr)
{
    forAll(shapeTypes_, i)
    {
        if (shapeTypes_[i] != "splitHex")
        {
            cellModels_[i] = cellModeller::lookup(shapeTypes_[i]);

            if (!cellModels_[i])
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown cell shape " << shapeTypes_[i]
                    << exit(FatalIOError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::cellShapes::~cellShapes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::cellShapes::generate() const
{
    boolList selectedCells(mesh_.nCells(), false);

    if (cellModels_.size() == 1 && !cellModels_[0])
    {
        for (label celli=0; celli<mesh_.nCells(); celli++)
        {
            const cellFeatures superCell(mesh_, featureCos, celli);

            if (hexMatcher().isA(superCell.faces()))
            {
                selectedCells[celli] = true;
            }
        }
    }
    else
    {
        const cellShapeList& cellShapes = mesh_.cellShapes();

        forAll(cellShapes, celli)
        {
            forAll(cellModels_, i)
            {
                if (cellModels_[i])
                {
                    if (cellShapes[celli].model() == *cellModels_[i])
                    {
                        selectedCells[celli] = true;
                        break;
                    }
                }
                else
                {
                    const cellFeatures superCell(mesh_, featureCos, celli);

                    if (hexMatcher().isA(superCell.faces()))
                    {
                        selectedCells[celli] = true;
                        break;
                    }
                }
            }
        }
    }

    return zoneSet
    (
        new cellZone
        (
            zoneName_,
            indices(selectedCells),
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
