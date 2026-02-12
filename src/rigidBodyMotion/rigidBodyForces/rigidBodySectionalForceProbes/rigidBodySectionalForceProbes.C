/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "rigidBodySectionalForceProbes.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rigidBodySectionalForceProbes, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        rigidBodySectionalForceProbes,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::functionObjects::rigidBodySectionalForceProbes::distances() const
{
    tmp<scalarField> tdistances =
        scalarField(coordinates_, coordinateOrders_)
      - localOrigin()[axisi()];
    scalarField& distances = tdistances.ref();

    const scalar maxDistance = gMax(patchPointDistances());

    distances.append((1 + 2*small)*max(maxDistance, distances.last()));

    return tdistances;
}


void Foam::functionObjects::rigidBodySectionalForceProbes::writeFileHeader
(
    const label
)
{
    forAll(coordinates_, i)
    {
        writeHeaderValue
        (
            file(),
            "Coordinate " + Foam::name(i),
            coordinates_[i]
        );
    }

    const Foam::Omanip<int> w = valueWidth(1);

    writeCommented(file(), "Coordinate");

    forAll(coordinates_, i)
    {
        file()
            << w << i << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' ';
    }
    file().endl();

    writeCommented(file(), "Time");

    forAll(coordinates_, i)
    {
        file()
            << w << "Fluid Force" << w << ' ' << w << ' '
            << w << "Body Force" << w << ' ' << w << ' '
            << w << "Total Force" << w << ' ' << w << ' '
            << w << "Fluid Moment" << w << ' ' << w << ' '
            << w << "Body Moment" << w << ' ' << w << ' '
            << w << "Total Moment" << w << ' ' << w << ' ';
    }
    file().endl();

    writeCommented(file(), "");

    forAll(coordinates_, i)
    {
        file()
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z';
    }
    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodySectionalForceProbes::
rigidBodySectionalForceProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    rigidBodySectionalForcesBase(name, runTime, dict),
    logFiles(mesh(), name),
    coordinates_(),
    coordinateOrders_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodySectionalForceProbes::
~rigidBodySectionalForceProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rigidBodySectionalForceProbes::read
(
    const dictionary& dict
)
{
    rigidBodySectionalForcesBase::read(dict);

    coordinates_ = dict.lookup<List<scalar>>("coordinates");

    if (coordinates_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "At least one coordinate must be specified"
            << exit(FatalIOError);
    }

    coordinateOrders_.setSize(coordinates_.size());
    sortedOrder(coordinates_, coordinateOrders_);

    resetName(typeName);

    return true;
}


bool Foam::functionObjects::rigidBodySectionalForceProbes::write()
{
    logFiles::write();

    vectorField fluidForce(coordinates_.size() + 1, Zero);
    vectorField bodyForce(coordinates_.size() + 1, Zero);
    vectorField totalForce(coordinates_.size() + 1, Zero);
    vectorField fluidMoment(coordinates_.size() + 1, Zero);
    vectorField bodyMoment(coordinates_.size() + 1, Zero);
    vectorField totalMoment(coordinates_.size() + 1, Zero);

    // Calculate the fluid contribution
    rigidBodySectionalForcesBase::addFluid(fluidForce, fluidMoment);

    // Calculate the body contribution
    rigidBodySectionalForcesBase::addBody(bodyForce, bodyMoment);

    // Sum to create total sectional forces and moments
    totalForce = fluidForce + bodyForce;
    totalMoment = fluidMoment + bodyMoment;

    // Reorder to match that of the input coordinates
    fluidForce.map(fluidForce, coordinateOrders_);
    bodyForce.map(bodyForce, coordinateOrders_);
    totalForce.map(totalForce, coordinateOrders_);
    fluidMoment.map(fluidMoment, coordinateOrders_);
    bodyMoment.map(bodyMoment, coordinateOrders_);
    totalMoment.map(totalMoment, coordinateOrders_);

    // Write out the graph
    if (Pstream::master())
    {
        writeTime(file());

        const Foam::Omanip<int> w = valueWidth(1);

        forAll(coordinates_, i)
        {
            file()
                << w << fluidForce[i].x()
                << w << fluidForce[i].y()
                << w << fluidForce[i].z()
                << w << bodyForce[i].x()
                << w << bodyForce[i].y()
                << w << bodyForce[i].z()
                << w << totalForce[i].x()
                << w << totalForce[i].y()
                << w << totalForce[i].z()
                << w << fluidMoment[i].x()
                << w << fluidMoment[i].y()
                << w << fluidMoment[i].z()
                << w << bodyMoment[i].x()
                << w << bodyMoment[i].y()
                << w << bodyMoment[i].z()
                << w << totalMoment[i].x()
                << w << totalMoment[i].y()
                << w << totalMoment[i].z();
        }

        file().endl();
    }

    return true;
}


// ************************************************************************* //
