/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fvMotionSolverEngineMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMotionSolverEngineMesh, 0);
    addToRunTimeSelectionTable(engineMesh, fvMotionSolverEngineMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolverEngineMesh::fvMotionSolverEngineMesh(const IOobject& io)
:
    engineMesh(io),
    pistonLayers_("pistonLayers", dimLength, 0.0),
    motionSolver_
    (
        *this,
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            engineDB_.engineDict()
        )
    )
{
    engineDB_.engineDict().readIfPresent("pistonLayers", pistonLayers_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolverEngineMesh::~fvMotionSolverEngineMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMotionSolverEngineMesh::move()
{
    scalar deltaZ = engineDB_.pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << endl;

    // Position of the top of the static mesh layers above the piston
    scalar pistonPlusLayers = pistonPosition_.value() + pistonLayers_.value();

    scalar pistonSpeed = deltaZ/engineDB_.deltaTValue();

    motionSolver_.pointMotionU().boundaryField()[pistonIndex_] == pistonSpeed;

    {
        scalarField linerPoints
        (
            boundary()[linerIndex_].patch().localPoints().component(vector::Z)
        );

        motionSolver_.pointMotionU().boundaryField()[linerIndex_] ==
            pistonSpeed*pos(deckHeight_.value() - linerPoints)
           *(deckHeight_.value() - linerPoints)
           /(deckHeight_.value() - pistonPlusLayers);
    }

    motionSolver_.solve();

    if (engineDB_.foundObject<surfaceScalarField>("phi"))
    {
        surfaceScalarField& phi =
            const_cast<surfaceScalarField&>
            (engineDB_.lookupObject<surfaceScalarField>("phi"));

        const volScalarField& rho =
            engineDB_.lookupObject<volScalarField>("rho");

        const volVectorField& U =
            engineDB_.lookupObject<volVectorField>("U");

        bool absolutePhi = false;
        if (moving())
        {
            phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
            absolutePhi = true;
        }

        movePoints(motionSolver_.curPoints());

        if (absolutePhi)
        {
            phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        }
    }
    else
    {
        movePoints(motionSolver_.curPoints());
    }


    pistonPosition_.value() += deltaZ;

    Info<< "clearance: " << deckHeight_.value() - pistonPosition_.value() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;
}


// ************************************************************************* //
