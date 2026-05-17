/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotion_pointMeshMover.H"
#include "timeIOdictionary.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pointMeshMovers
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotion, 0);

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        sixDoFRigidBodyMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::septernion>
Foam::pointMeshMovers::sixDoFRigidBodyMotion::transforms0() const
{
    // Assume the external body is stationary
    return List<septernion>
    (
        {sixDoFRigidBodyMotion::transform0(), septernion::I}
    );
}


void Foam::pointMeshMovers::sixDoFRigidBodyMotion::moveBodies()
{
    const Time& t = poly().time();

    if (poly().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << poly().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-stepbool
    bool firstIter = false;
    if (curTimeIndex_ != t.timeIndex())
    {
        newTime();
        curTimeIndex_ = t.timeIndex();
        firstIter = true;
    }

    dimensionedVector g(g_);

    if (poly().foundObject<uniformDimensionedVectorField>("g"))
    {
        g = poly().lookupObject<uniformDimensionedVectorField>("g");
    }

    // scalar ramp = min(max((t.value() - 5)/10, 0), 1);
    scalar ramp = 1.0;

    if (test_)
    {
        update
        (
            firstIter,
            ramp*(mass()*g.value()),
            ramp*(mass()*(momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
    else
    {
        functionObjects::forces f
        (
            functionObjects::forces::typeName,
            t,
            dictionary::entries
            (
                "type", functionObjects::forces::typeName,
                "patches", bodyMeshes_[0].patches(),
                "rhoInf", rhoInf_,
                "rho", rhoName_,
                "CofR", centreOfRotation()
            )
        );

        f.calcForcesMoments();

        update
        (
            firstIter,
            ramp*(f.forceEff() + mass()*g.value()),
            ramp
           *(
               f.momentEff()
             + mass()*(momentArm() ^ g.value())
            ),
            t.deltaTValue(),
            t.deltaT0Value()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMeshMovers::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    pointMeshMovers::multiRigidBody(mesh, dict),
    Foam::sixDoFRigidBodyMotion
    (
        dict,
        typeIOobject<timeIOdictionary>
        (
            "sixDoFRigidBodyMotionState",
            mesh.time().name(),
            "uniform",
            mesh
        ).headerOk()
      ? timeIOdictionary
        (
            IOobject
            (
                "sixDoFRigidBodyMotionState",
                mesh.time().name(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : dict
    ),
    test_(dict.lookupOrDefault<Switch>("test", false)),
    rhoInf_(1.0),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    g_("g", dimAcceleration, dict, vector::zero),
    curTimeIndex_(-1)
{
    if (rhoName_ == "rhoInf")
    {
        rhoInf_ = dict.lookup<scalar>("rhoInf");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointMeshMovers::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointMeshMovers::sixDoFRigidBodyMotion::write() const
{
    timeIOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            poly().time().name(),
            "uniform",
            poly(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    state().write(dict);

    return
        pointMeshMovers::multiRigidBody::write()
     && dict.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            poly().time().writeCompression(),
            true
        );
}


// ************************************************************************* //
