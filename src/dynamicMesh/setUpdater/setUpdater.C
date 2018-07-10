/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "setUpdater.H"
#include "polyTopoChanger.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setUpdater, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        setUpdater,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::setUpdater::setUpdater
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, Switch(dict.lookup("active")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setUpdater::~setUpdater()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setUpdater::changeTopology() const
{
    // I am never cause of changeTopo
    return false;
}


void Foam::setUpdater::setRefinement(polyTopoChange&) const
{}


void Foam::setUpdater::modifyMotionPoints(pointField&) const
{}


void Foam::setUpdater::updateMesh(const mapPolyMesh& morphMap)
{
    // Mesh has changed topologically. Update all sets.
    if (debug)
    {
        Pout<< "setUpdater::updateMesh(const mapPolyMesh& morphMap)"
            << endl;
    }

    updateSets<cellSet>(morphMap);
    updateSets<faceSet>(morphMap);
    updateSets<pointSet>(morphMap);
}


void Foam::setUpdater::write(Ostream& os) const
{
    os  << nl << type() << nl;
}


void Foam::setUpdater::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
