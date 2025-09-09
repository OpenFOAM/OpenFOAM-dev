/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "periodic.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(periodic, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            periodic,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::zoneGenerators::periodic::activate() const
{
    const scalar t =
        repeat_ > 0
      ? begin_ + fmod(mesh_.time().userTimeValue() - begin_, repeat_)
      : mesh_.time().userTimeValue();

    const bool active = t >= begin_ && t < end_;

    return deactivate_ ? !active : active;
}


bool Foam::zoneGenerators::periodic::stateChanged() const
{
    if (activate())
    {
        return !active_;
    }
    else
    {
        return active_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::periodic::periodic
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerator_(zoneGenerator::New(mesh, dict)),
    begin_
    (
        dict.lookupOrDefault<scalar>
        (
            "begin",
            unitNone,
            mesh.time().beginTime().value()
        )

    ),
    end_
    (
        dict.lookupOrDefault<scalar>
        (
            "end",
            unitNone,
            mesh.time().endTime().value()
        )

    ),
    repeat_(dict.lookupOrDefault<scalar>("repeat", unitNone, 0)),
    deactivate_(dict.lookupOrDefault<Switch>("deactivate", false)),
    active_(false)
{
    moveUpdate_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::periodic::~periodic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::periodic::generate() const
{
    if (activate())
    {
        active_ = true;
        return zoneGenerator_->generate().clone(zoneName_);
    }
    else
    {
        active_ = false;

        if (!nullZoneSet_.valid())
        {
            zoneSet zs(zoneGenerator_->generate());

            nullZoneSet_ = zoneSet
            (
                zs.pValid()
              ? new pointZone
                (
                    zoneName_,
                    labelList::null(),
                    mesh_.pointZones(),
                    moveUpdate_,
                    true
                )
              : nullptr,
                zs.cValid()
              ? new cellZone
                (
                    zoneName_,
                    labelList::null(),
                    mesh_.cellZones(),
                    moveUpdate_,
                    true
                )
              : nullptr,
                zs.fValid()
              ? new faceZone
                (
                    zoneName_,
                    labelList::null(),
                    boolList::null(),
                    mesh_.faceZones(),
                    moveUpdate_,
                    true
                )
              : nullptr
            );
        }

        // Return a zoneSet of references to nullZoneSet_
        return zoneSet(nullZoneSet_, false);
    }
}


Foam::zoneSet Foam::zoneGenerators::periodic::movePoints() const
{
    // Regenerate the zones if they update on move or the state has changed
    if (zoneGenerator_->moveUpdate() || stateChanged())
    {
        return generate();
    }
    else
    {
        return zoneSet();
    }
}


// ************************************************************************* //
