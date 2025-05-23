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

Class
    Foam::zoneGenerators::periodic

Description
    A zoneGenerator which activates the zoneSet returned by the given
    zoneGenerator for a given period with optional repetition.

Usage
    \table
        Property  | Description                   | Required | Default value
        type      | Type: periodic                | yes      |
        begin     | Begin time for the zones      | no       | Time::beginTime()
        end       | End time for the zones        | no       | Time::endTime()
        repeat    | Repetition period             | no       |
    \endtable
    These options are followed by a single zoneGenerator specification.

    Example of the specification for a periodic box cellZone in the
    constant/zonesGenerator file:
    \verbatim
        box1
        {
            type    periodic;

            begin   0.001;
            end     0.002;
            repeat  0.002;

            box1
            {
                type        box;
                zoneType    cell;
                moveUpdate  true;

                min     (-0.0075 0 -1);
                max     (-0.003 0.0025 1);
            }
        }
    \endverbatim

SourceFiles
    periodic.C

\*---------------------------------------------------------------------------*/

#ifndef periodic_H
#define periodic_H

#include "zoneGeneratorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace zoneGenerators
{

/*---------------------------------------------------------------------------*\
                          Class periodic Declaration
\*---------------------------------------------------------------------------*/

class periodic
:
    public zoneGenerator
{
    // Private Data

        autoPtr<zoneGenerator> zoneGenerator_;

        //- Optional begin time for the zones being active
        //  Defaults to Time::beginTime
        const scalar begin_;

        //- Optional end time for the zones being active
        //  Defaults to Time::endTime
        const scalar end_;

        //- Optional repetition period for the zones being active
        //  Defaults to no repetition
        const scalar repeat_;

        //- Active state
        mutable bool active_;

        //- Cached null zoneSet
        mutable zoneSet nullZoneSet_;


    // Private Member Functions

        //- Returns true when the zones should be activated
        bool activate() const;

        //- Returns true if the active state of the zones has changed
        bool stateChanged() const;


public:

    //- Runtime type information
    TypeName("periodic");


    // Constructors

        //- Construct from dictionary
        periodic
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~periodic();


    // Member Functions

        virtual zoneSet generate() const;

        virtual zoneSet movePoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace zoneGenerators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
