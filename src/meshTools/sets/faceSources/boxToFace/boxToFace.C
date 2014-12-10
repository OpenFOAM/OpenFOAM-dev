/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "boxToFace.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxToFace, 0);

addToRunTimeSelectionTable(topoSetSource, boxToFace, word);

addToRunTimeSelectionTable(topoSetSource, boxToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::boxToFace::usage_
(
    boxToFace::typeName,
    "\n    Usage: boxToFace ((minx miny minz) (maxx maxy maxz))\n\n"
    "    Select all face with faceCentre within bounding box\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boxToFace::combine(topoSet& set, const bool add) const
{
    const pointField& ctrs = mesh_.faceCentres();

    forAll(ctrs, faceI)
    {
        forAll(bbs_, i)
        {
            if (bbs_[i].contains(ctrs[faceI]))
            {
                addOrDelete(set, faceI, add);
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    const treeBoundBoxList& bbs
)
:
    topoSetSource(mesh),
    bbs_(bbs)
{}


// Construct from dictionary
Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    bbs_
    (
        dict.found("box")
      ? treeBoundBoxList(1, treeBoundBox(dict.lookup("box")))
      : dict.lookup("boxes")
    )
{}


// Construct from Istream
Foam::boxToFace::boxToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    bbs_(1, treeBoundBox(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boxToFace::~boxToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boxToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces with centre within boxes " << bbs_ << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces with centre within boxes " << bbs_ << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
