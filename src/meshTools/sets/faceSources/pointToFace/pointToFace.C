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

#include "pointToFace.H"
#include "polyMesh.H"
#include "pointSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, pointToFace, word);
    addToRunTimeSelectionTable(topoSetSource, pointToFace, istream);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::pointToFace::pointAction,
        3
    >::names[] =
    {
        "any",
        "all",
        "edge"
    };
}


Foam::topoSetSource::addToUsageTable Foam::pointToFace::usage_
(
    pointToFace::typeName,
    "\n    Usage: pointToFace <pointSet> any|all|edge\n\n"
    "    Select faces with\n"
    "    -any point in the pointSet\n"
    "    -all points in the pointSet\n\n"
    "    -two consecutive points (an edge) in the pointSet\n\n"
);

const Foam::NamedEnum<Foam::pointToFace::pointAction, 3>
    Foam::pointToFace::pointActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointToFace::combine(topoSet& set, const bool add) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName_);

    if (option_ == ANY)
    {
        // Add faces with any point in loadedSet
        forAllConstIter(pointSet, loadedSet, iter)
        {
            const label pointI = iter.key();
            const labelList& pFaces = mesh_.pointFaces()[pointI];

            forAll(pFaces, pFaceI)
            {
                addOrDelete(set, pFaces[pFaceI], add);
            }
        }
    }
    else if (option_ == ALL)
    {
        // Add all faces whose points are all in set.

        // Count number of points using face.
        Map<label> numPoints(loadedSet.size());

        forAllConstIter(pointSet, loadedSet, iter)
        {
            const label pointI = iter.key();
            const labelList& pFaces = mesh_.pointFaces()[pointI];

            forAll(pFaces, pFaceI)
            {
                const label faceI = pFaces[pFaceI];

                Map<label>::iterator fndFace = numPoints.find(faceI);

                if (fndFace == numPoints.end())
                {
                    numPoints.insert(faceI, 1);
                }
                else
                {
                    fndFace()++;
                }
            }
        }


        // Include faces that are referenced as many times as there are points
        // in face -> all points of face
        forAllConstIter(Map<label>, numPoints, iter)
        {
            const label faceI = iter.key();

            if (iter() == mesh_.faces()[faceI].size())
            {
                addOrDelete(set, faceI, add);
            }
        }
    }
    else if (option_ == EDGE)
    {
        const faceList& faces = mesh_.faces();
        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, fp)
            {
                if (loadedSet.found(f[fp]) && loadedSet.found(f.nextLabel(fp)))
                {
                    addOrDelete(set, faceI, add);
                    break;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(pointActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::pointToFace::pointToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointToFace::~pointToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding faces according to pointSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing faces according to pointSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
