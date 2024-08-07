/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "rigidBodyPoints.H"
#include "motionSolver_fvMeshMover.H"
#include "motionSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rigidBodyPoints, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        rigidBodyPoints,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::RBD::rigidBodyMotion&
Foam::functionObjects::rigidBodyPoints::motion() const
{
    const fvMeshMovers::motionSolver& mover =
        refCast<const fvMeshMovers::motionSolver>(mesh_.mover());

    return (refCast<const RBD::rigidBodyMotion>(mover.motion()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodyPoints::rigidBodyPoints
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    angularVelocityUnits_("[rad/s]"),
    angularAccelerationUnits_("[rad/s^2]")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodyPoints::~rigidBodyPoints()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rigidBodyPoints::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    angularVelocityUnits_.readIfPresent("angularVelocityUnits", dict);
    angularAccelerationUnits_.readIfPresent("angularAccelerationUnits", dict);

    dict.lookup("body") >> body_;

    const HashTable<point> pointsTable(dict.lookup("points"));
    names_.setSize(pointsTable.size());
    points_.setSize(pointsTable.size());

    label i = 0;
    forAllConstIter(HashTable<point>, pointsTable, iter)
    {
        names_[i] = iter.key();
        points_[i++] = iter();
    }

    resetNames(names_);

    return true;
}


void Foam::functionObjects::rigidBodyPoints::writeFileHeader(const label i)
{
    OFstream& file = this->files()[i];

    writeHeader(file, "Body point motion");
    writeHeaderValue(file, "Body", body_);
    writeHeaderValue(file, "Point", points_[i]);
    writeHeaderValue(file, "Angular Velocity Units", angularVelocityUnits_);
    writeHeaderValue
    (
        file,
        "Angular Acceleration Units",
        angularAccelerationUnits_
    );
    writeCommented(file, "Time");

    file<< tab
        << "Position" << tab
        << "Linear velocity" << tab
        << "Angular velocity" << tab
        << "Linear acceleration" << tab
        << "Angular acceleration" << endl;
}


bool Foam::functionObjects::rigidBodyPoints::execute()
{
    return true;
}


bool Foam::functionObjects::rigidBodyPoints::write()
{
    logFiles::write();

    if (Pstream::master())
    {
        const RBD::rigidBodyMotion& motion = this->motion();

        const label bodyID = motion.bodyIndex(body_);

        forAll(points_, i)
        {
            const vector p(motion.p(bodyID, points_[i]));
            const spatialVector v(motion.v(bodyID, points_[i]));
            const spatialVector a(motion.a(bodyID, points_[i]));

            writeTime(files()[i]);
            files()[i]
                << tab
                << p  << tab
                << v.l() << tab
                << angularVelocityUnits_.toUser(v.w()) << tab
                << a.l() << tab
                << angularAccelerationUnits_.toUser(a.w()) << endl;
        }
    }

    return true;
}


// ************************************************************************* //
