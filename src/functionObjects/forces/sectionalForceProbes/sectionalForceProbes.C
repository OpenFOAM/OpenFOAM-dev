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

#include "sectionalForceProbes.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sectionalForceProbes, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        sectionalForceProbes,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::vector Foam::functionObjects::sectionalForceProbes::normal() const
{
    return normal_;
}


Foam::point Foam::functionObjects::sectionalForceProbes::origin() const
{
    return points_.first();
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sectionalForceProbes::distances() const
{
    const scalar dot0 = points_.first() & normal_;
    const scalar dot1 = points_.last() & normal_;

    tmp<scalarField> tdistances =
        (pointField(points_, pointOrders_) & normal_) - dot0;
    scalarField& distances = tdistances.ref();

    const scalar maxDistance = gMax(patchPointDistances());

    distances.append((1 + 2*small)*max(maxDistance, dot1 - dot0));

    return tdistances;
}


void Foam::functionObjects::sectionalForceProbes::writeFileHeader(const label)
{
    forAll(points_, i)
    {
        writeHeaderValue
        (
            file(),
            "Point " + Foam::name(i),
            points_[i]
        );
    }

    const Foam::Omanip<int> w = valueWidth(1);

    writeCommented(file(), "Point");

    forAll(points_, i)
    {
        file()
            << w << i << w << ' ' << w << ' '
            << w << ' ' << w << ' ' << w << ' ';
    }
    file().endl();

    writeCommented(file(), "Time");

    forAll(points_, i)
    {
        file()
            << w << "Force" << w << ' ' << w << ' '
            << w << "Moment" << w << ' ' << w << ' ';
    }
    file().endl();

    writeCommented(file(), "");

    forAll(points_, i)
    {
        file()
            << w << 'x' << w << 'y' << w << 'z'
            << w << 'x' << w << 'y' << w << 'z';
    }
    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForceProbes::sectionalForceProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sectionalForcesBase(name, runTime, dict),
    logFiles(mesh(), name),
    points_(),
    pointOrders_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForceProbes::~sectionalForceProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sectionalForceProbes::read(const dictionary& dict)
{
    sectionalForcesBase::read(dict);

    normal_ = normalised(dict.lookup<vector>("normal"));
    points_ = dict.lookup<List<point>>("points");

    if (points_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "At least one point must be specified"
            << exit(FatalIOError);
    }

    // Check the points are in a line
    SubField<vector> pointsNot0(points_, points_.size() - 1, 1);
    const scalarField pointDeviation
    (
        mag((tensor::I - sqr(normal_)) & (pointsNot0 - points_.first()))
       /mag(points_.last() - points_.first())
    );
    const scalar maxPointDeviation = max(pointDeviation);
    if (maxPointDeviation > rootSmall)
    {
        WarningInFunction
            << "Points of " << typeName << " function " << name()
            << " are out of line in the direction " << normal_ << " by up to "
            << maxPointDeviation*100 << "\% of the total span. This means"
            << " that the moments computed at the points have different"
            << " origins, and therefore may not be comparable." << endl;
    }

    pointOrders_.setSize(points_.size());
    sortedOrder((points_ & normal_)(), pointOrders_);

    resetName(typeName);

    return true;
}


bool Foam::functionObjects::sectionalForceProbes::write()
{
    logFiles::write();

    vectorField force(points_.size() + 1, Zero);
    vectorField moment(points_.size() + 1, Zero);

    sectionalForcesBase::addFluid(force, moment);

    force.map(force, pointOrders_);
    moment.map(moment, pointOrders_);

    // Write out the graph
    if (Pstream::master())
    {
        writeTime(file());

        const Foam::Omanip<int> w = valueWidth(1);

        forAll(points_, i)
        {
            file()
                << w << force[i].x()
                << w << force[i].y()
                << w << force[i].z()
                << w << moment[i].x()
                << w << moment[i].y()
                << w << moment[i].z();
        }

        file().endl();
    }

    return true;
}


// ************************************************************************* //
