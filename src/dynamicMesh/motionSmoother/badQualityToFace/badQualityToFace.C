/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "badQualityToFace.H"
#include "polyMesh.H"
#include "motionSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(badQualityToFace, 0);

addToRunTimeSelectionTable(topoSetSource, badQualityToFace, word);

addToRunTimeSelectionTable(topoSetSource, badQualityToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::badQualityToFace::usage_
(
    badQualityToFace::typeName,
    "\n    Usage: badQualityToFace mesh-quality-dictionary\n\n"
    "    Select all faces that do not satisfy the selection criterion\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::badQualityToFace::combine(topoSet& set, const bool add) const
{
    faceSet faces(mesh_, "meshQualityFaces", mesh_.nFaces()/100+1);
    motionSmoother::checkMesh(false, mesh_, dict_, faces);
    faces.sync(mesh_);

    forAllConstIter(faceSet, faces, iter)
    {
        label faceI = iter.key();
        addOrDelete(set, faceI, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::badQualityToFace::badQualityToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    dict_(dict)
{}


// Construct from Istream
Foam::badQualityToFace::badQualityToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    dict_(is)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::badQualityToFace::~badQualityToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::badQualityToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding bad-quality faces" << endl;
        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing bad-quality faces" << endl;
        combine(set, false);
    }
}


// ************************************************************************* //
