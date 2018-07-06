/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "externalDisplacementMeshMover.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalDisplacementMeshMover, 0);
    defineRunTimeSelectionTable(externalDisplacementMeshMover, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover::externalDisplacementMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement
)
:
    baffles_(baffles),
    pointDisplacement_(pointDisplacement)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::externalDisplacementMeshMover>
Foam::externalDisplacementMeshMover::New
(
    const word& type,
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement
)
{
    Info<< "Selecting externalDisplacementMeshMover " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown externalDisplacementMeshMover type "
            << type << nl << nl
            << "Valid externalDisplacementMeshMover types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<externalDisplacementMeshMover>
    (
        cstrIter()(dict, baffles, pointDisplacement)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover::~externalDisplacementMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalDisplacementMeshMover::movePoints(const pointField&)
{
    // No local data to update
}


void Foam::externalDisplacementMeshMover::updateMesh(const mapPolyMesh& mpm)
{
    // Renumber baffles
    DynamicList<labelPair> newBaffles(baffles_.size());
    forAll(baffles_, i)
    {
        label f0 = mpm.reverseFaceMap()[baffles_[i].first()];
        label f1 = mpm.reverseFaceMap()[baffles_[i].second()];

        if (f0 >= 0 && f1 >= 0)
        {
            newBaffles.append(labelPair(f0, f1));
        }
    }
    newBaffles.shrink();
    baffles_.transfer(newBaffles);
}


// ************************************************************************* //
