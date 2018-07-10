/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "sixDoFRigidBodyState.H"
#include "dynamicMotionSolverFvMesh.H"
#include "sixDoFRigidBodyMotionSolver.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sixDoFRigidBodyState, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sixDoFRigidBodyState,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sixDoFRigidBodyState::sixDoFRigidBodyState
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sixDoFRigidBodyState::~sixDoFRigidBodyState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sixDoFRigidBodyState::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    angleFormat_ = dict.lookupOrDefault<word>("angleFormat", "radians");

    return true;
}


void Foam::functionObjects::sixDoFRigidBodyState::writeFileHeader(const label)
{
    OFstream& file = this->file();

    writeHeader(file, "Motion State");
    writeHeaderValue(file, "Angle Units", angleFormat_);
    writeCommented(file, "Time");

    file<< tab
        << "centreOfRotation" << tab
        << "centreOfMass" << tab
        << "rotation" << tab
        << "velocity" << tab
        << "omega" << endl;
}


bool Foam::functionObjects::sixDoFRigidBodyState::execute()
{
    return true;
}


bool Foam::functionObjects::sixDoFRigidBodyState::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        const dynamicMotionSolverFvMesh& mesh =
            refCast<const dynamicMotionSolverFvMesh>(obr_);

        const sixDoFRigidBodyMotionSolver& motionSolver_ =
            refCast<const sixDoFRigidBodyMotionSolver>(mesh.motion());

        const sixDoFRigidBodyMotion& motion = motionSolver_.motion();

        vector rotationAngle
        (
            quaternion(motion.orientation()).eulerAngles(quaternion::XYZ)
        );

        vector angularVelocity(motion.omega());

        if (angleFormat_ == "degrees")
        {
            rotationAngle.x() = radToDeg(rotationAngle.x());
            rotationAngle.y() = radToDeg(rotationAngle.y());
            rotationAngle.z() = radToDeg(rotationAngle.z());

            angularVelocity.x() = radToDeg(angularVelocity.x());
            angularVelocity.y() = radToDeg(angularVelocity.y());
            angularVelocity.z() = radToDeg(angularVelocity.z());
        }

        writeTime(file());
        file()
            << tab
            << motion.centreOfRotation()  << tab
            << motion.centreOfMass()  << tab
            << rotationAngle  << tab
            << motion.v()  << tab
            << angularVelocity << endl;
    }

    return true;
}


// ************************************************************************* //
