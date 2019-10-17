/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::crankConRod

Description
    Manage time in terms of engine RPM and crank-angle.

    When crankConRod is in effect, the userTime is reported in degrees
    crank-angle instead of in seconds. The RPM to be used is specified in
    \c constant/engineGeometry. If only a time conversion is required,
    the geometric engine parameters can be dropped or set to zero.

    For example,
    \verbatim
        rpm             rpm  [0 0 -1 0 0]  2000;

        conRodLength    conRodLength  [0 1 0 0 0] 0.0;
        bore            bore          [0 1 0 0 0] 0.0;
        stroke          stroke        [0 1 0 0 0] 0.0;
        clearance       clearance     [0 1 0 0 0] 0.0;
    \endverbatim

SourceFiles
    crankConRod.C

\*---------------------------------------------------------------------------*/

#ifndef crankConRod_H
#define crankConRod_H

#include "engineTime.H"
#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class crankConRod Declaration
\*---------------------------------------------------------------------------*/

class crankConRod
:
    public engineTime
{
    // Private Data

        //- RPM is required
        dimensionedScalar rpm_;

        //- Optional engine geometry parameters
        dimensionedScalar conRodLength_;
        dimensionedScalar bore_;
        dimensionedScalar stroke_;
        dimensionedScalar clearance_;


    // Private Member Functions

        //- Adjust read time values
        void timeAdjustment();


public:

    //- Runtime type information
    TypeName("crankConRod");


    // Constructors

        //- Construct from objectRegistry arguments
        crankConRod
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName = "system",
            const fileName& constantName = "constant",
            const fileName& dictName = "engineGeometry"
        );

        //- Disallow default bitwise copy construction
        crankConRod(const crankConRod&) = delete;


    //- Destructor
    virtual ~crankConRod()
    {}


    // Member Functions

        // Conversion

            //- Convert degrees to seconds (for given engine speed in RPM)
            scalar degToTime(const scalar theta) const;

            //- Convert seconds to degrees (for given engine speed in RPM)
            scalar timeToDeg(const scalar t) const;

            //- Calculate the piston position from the engine geometry
            //  and given crank angle.
            scalar pistonPosition(const scalar theta) const;


        // Access

            //- Return the engines current operating RPM
            const dimensionedScalar& rpm() const
            {
                return rpm_;
            }

            //- Return the engines connecting-rod length
            const dimensionedScalar& conRodLength() const
            {
                return conRodLength_;
            }

            //- Return the engines bore
            const dimensionedScalar& bore() const
            {
                return bore_;
            }

            //- Return the engines stroke
            const dimensionedScalar& stroke() const
            {
                return stroke_;
            }

            //- Return the engines clearance-gap
            const dimensionedScalar& clearance() const
            {
                return clearance_;
            }


            //- Return current crank-angle
            virtual scalar theta() const;

            //- Return time unit
            virtual word unit() const;

            //- Return current crank-angle translated to a single revolution
            //  (value between -180 and 180 with 0 = top dead centre)
            scalar thetaRevolution() const;

            //- Return crank-angle increment
            virtual scalar deltaTheta() const;


        // Member Functions overriding the virtual functions in time

            //- Convert the user-time (CA deg) to real-time (s).
            virtual scalar userTimeToTime(const scalar theta) const;

            //- Convert the real-time (s) into user-time (CA deg)
            virtual scalar timeToUserTime(const scalar t) const;

            //- Ratio between real-time and user-time
            virtual scalar timeToUserTimeRatio() const;

            //- Read the control dictionary and set the write controls etc.
            virtual void readDict();


        // Edit

            //- Read the controlDict and set all the parameters
            virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const crankConRod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
