/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "rigidBodyMotion.H"
#include "rigidBodySolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyMotion::rigidBodyMotion()
:
    rigidBodyModel(),
    motionState_(*this),
    motionState0_(*this),
    aRelax_(1.0),
    aDamp_(1.0),
    report_(false),
    solver_(NULL)
{}


Foam::RBD::rigidBodyMotion::rigidBodyMotion
(
    const dictionary& dict
)
:
    rigidBodyModel(dict),
    motionState_(*this),
    motionState0_(motionState_),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    aDamp_(dict.lookupOrDefault<scalar>("accelerationDamping", 1.0)),
    report_(dict.lookupOrDefault<Switch>("report", false)),
    solver_(rigidBodySolver::New(*this, dict.subDict("solver")))
{}


Foam::RBD::rigidBodyMotion::rigidBodyMotion
(
    const dictionary& dict,
    const dictionary& stateDict
)
:
    rigidBodyModel(dict),
    motionState_(*this, stateDict),
    motionState0_(motionState_),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    aDamp_(dict.lookupOrDefault<scalar>("accelerationDamping", 1.0)),
    report_(dict.lookupOrDefault<Switch>("report", false)),
    solver_(rigidBodySolver::New(*this, dict.subDict("solver")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::rigidBodyMotion::~rigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::rigidBodyMotion::update
(
    scalar deltaT,
    scalar deltaT0,
    const scalarField& tau,
    const Field<spatialVector>& fx
)
{
    if (Pstream::master())
    {
        solver_->solve(deltaT, deltaT0, tau, fx);

        if (report_)
        {
            status();
        }
    }

    Pstream::scatter(motionState_);
}


void Foam::RBD::rigidBodyMotion::status() const
{
    /*
    Info<< "Rigid body motion" << nl
        << "    Centre of rotation: " << centreOfRotation() << nl
        << "    Centre of mass: " << centreOfMass() << nl
        << "    Orientation: " << orientation() << nl
        << "    Linear velocity: " << v() << nl
        << "    Angular velocity: " << omega()
        << endl;
    */
}

/*
Foam::tmp<Foam::pointField> Foam::RBD::rigidBodyMotion::transform
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfRotation_))
    );
}


Foam::tmp<Foam::pointField> Foam::RBD::rigidBodyMotion::transform
(
    const pointField& initialPoints,
    const scalarField& scale
) const
{
    // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfRotation() - initialCentreOfRotation(),
        quaternion(Q() & initialQ().T())
    );

    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints.ref();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale[pointi] > 1 - SMALL)
            {
                points[pointi] = transform(initialPoints[pointi]);
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale[pointi]));

                points[pointi] =
                    initialCentreOfRotation()
                  + ss.transform
                    (
                        initialPoints[pointi]
                      - initialCentreOfRotation()
                    );
            }
        }
    }

    return tpoints;
}
*/

// ************************************************************************* //
