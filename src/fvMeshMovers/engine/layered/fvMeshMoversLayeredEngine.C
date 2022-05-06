/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshMoversLayeredEngine.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(layeredEngine, 0);
    addToRunTimeSelectionTable(fvMeshMover, layeredEngine, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::layeredEngine::layeredEngine(fvMesh& mesh)
:
    engine(mesh),
    pistonLayers_("pistonLayers", dimLength, meshCoeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::layeredEngine::~layeredEngine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshMovers::layeredEngine::update()
{
    const scalar deltaZ = pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;

    // Position of the top of the static mesh layers above the piston
    const scalar pistonPlusLayers =
        pistonPosition_.value() + pistonLayers_.value();

    pointField newPoints(mesh().points());

    forAll(newPoints, pointi)
    {
        point& p = newPoints[pointi];

        if (p.z() < pistonPlusLayers)           // In piston bowl
        {
            p.z() += deltaZ;
        }
        else if (p.z() < deckHeight_.value())   // In liner region
        {
            p.z() +=
                deltaZ
               *(deckHeight_.value() - p.z())
               /(deckHeight_.value() - pistonPlusLayers);
        }
    }

    mesh().movePoints(newPoints);

    pistonPosition_.value() += deltaZ;
    const scalar pistonSpeed = deltaZ/mesh().time().deltaTValue();

    Info<< "clearance: " << deckHeight_.value() - pistonPosition_.value() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;

    return true;
}


void Foam::fvMeshMovers::layeredEngine::topoChange(const polyTopoChangeMap&)
{
    NotImplemented;
}


void Foam::fvMeshMovers::layeredEngine::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshMovers::layeredEngine::distribute
(
    const polyDistributionMap&
)
{
    NotImplemented;
}


// ************************************************************************* //
