/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "rotatedBoxToCell.H"
#include "polyMesh.H"
#include "cellModeller.H"
#include "transform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotatedBoxToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, rotatedBoxToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rotatedBoxToCell::combine(topoSet& set, const bool add) const
{
    // Define a cell for the box
    const pointField boxPoints
    (
        {
            origin_,
            origin_ + i_,
            origin_ + i_ + j_,
            origin_ + j_,
            origin_ + k_,
            origin_ + k_ + i_,
            origin_ + k_ + i_ + j_,
            origin_ + k_ + j_
        }
    );

    const labelList boxVerts({0, 1, 2, 3, 4, 5, 6, 7});

    const cellModel& hex = *(cellModeller::lookup("hex"));

    // Get outwards pointing faces.
    const faceList boxFaces(cellShape(hex, boxVerts).faces());

    // Precalculate normals
    vectorField boxFaceNormals(boxFaces.size());
    forAll(boxFaces, i)
    {
        boxFaceNormals[i] = boxFaces[i].area(boxPoints);
    }

    // Check whether cell centre is inside all faces of box.

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        bool inside = true;

        forAll(boxFaces, i)
        {
            if
            (
                ((ctrs[celli] - boxPoints[boxFaces[i][0]]) & boxFaceNormals[i])
              > 0
            )
            {
                inside = false;
                break;
            }
        }

        if (inside)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatedBoxToCell::rotatedBoxToCell
(
    const polyMesh& mesh,
    const vector& origin,
    const vector& i,
    const vector& j,
    const vector& k
)
:
    topoSetSource(mesh),
    origin_(origin),
    i_(i),
    j_(j),
    k_(k)
{}


Foam::rotatedBoxToCell::rotatedBoxToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    origin_(),
    i_(),
    j_(),
    k_()
{
    if (dict.found("box"))
    {
        const boundBox bb(dict.lookup("box"));
        const vector c
        (
            dict.lookupOrDefault<point>("centre", dimLength, bb.midpoint())
        );
        const vector n1(normalised(dict.lookup<vector>("n1", dimless)));
        const vector n2(normalised(dict.lookup<vector>("n2", dimless)));

        const tensor R(rotationTensor(n1, n2));
        const pointField bbPoints(bb.points());

        origin_ = (R & (bb.min() - c)) + c;
        i_ = R & (bbPoints[1] - bb.min());
        j_ = R & (bbPoints[3] - bb.min());
        k_ = R & (bbPoints[4] - bb.min());
    }
    else
    {
        origin_ = dict.lookup<point>("origin", dimLength);
        i_ = dict.lookup<vector>("i", dimLength);
        j_ = dict.lookup<vector>("j", dimLength);
        k_ = dict.lookup<vector>("k", dimLength);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rotatedBoxToCell::~rotatedBoxToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatedBoxToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with center within rotated box " << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with center within rotated box " << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
