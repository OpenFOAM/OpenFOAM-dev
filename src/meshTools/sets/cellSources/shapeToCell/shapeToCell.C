/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "shapeToCell.H"
#include "polyMesh.H"
#include "unitConversion.H"
#include "hexMatcher.H"
#include "cellFeatures.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shapeToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, shapeToCell, word);
}


// Angle for polys to be considered splitHexes.
Foam::scalar Foam::shapeToCell::featureCos = Foam::cos(degToRad(10.0));


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shapeToCell::combine(topoSet& set, const bool add) const
{
    if (type_ == "splitHex")
    {
        for (label celli = 0; celli < mesh_.nCells(); celli++)
        {
            cellFeatures superCell(mesh_, featureCos, celli);

            if (hexMatcher().isA(superCell.faces()))
            {
                addOrDelete(set, celli, add);
            }
        }
    }
    else
    {
        const cellModel& wantedModel = *(cellModeller::lookup(type_));

        const cellShapeList& cellShapes = mesh_.cellShapes();

        forAll(cellShapes, celli)
        {
            if (cellShapes[celli].model() == wantedModel)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shapeToCell::shapeToCell
(
    const polyMesh& mesh,
    const word& type
)
:
    topoSetSource(mesh),
    type_(type)
{
    if (!cellModeller::lookup(type_) && (type_ != "splitHex"))
    {
        FatalErrorInFunction
            << "Illegal cell type " << type_ << exit(FatalError);
    }
}


Foam::shapeToCell::shapeToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    type_(dict.lookup("type"))
{
    if (!cellModeller::lookup(type_) && (type_ != "splitHex"))
    {
        FatalErrorInFunction
            << "Illegal cell type " << type_ << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shapeToCell::~shapeToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shapeToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells of type " << type_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of type " << type_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
