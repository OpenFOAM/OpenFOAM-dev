/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::engineTime

Description
    Manage time in terms of engine RPM and crank-angle.

    When engineTime is in effect, the userTime is reported in degrees
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

Note
   The engineTime can currently only be selected at compile-time.

SourceFiles
    engineTime.C

\*---------------------------------------------------------------------------*/

#ifndef engineTime_H
#define engineTime_H

#include "Time.H"
#include "dictionary.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class engineTime Declaration
\*---------------------------------------------------------------------------*/

class engineTime
:
    public Time
{
    // Private data

        IOdictionary dict_;

        //- RPM is required
        dimensionedScalar rpm_;

        //- Optional engine geometry parameters
        dimensionedScalar conRodLength_;
        dimensionedScalar bore_;
        dimensionedScalar stroke_;
        dimensionedScalar clearance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        engineTime(const engineTime&);

        //- Disallow default bitwise assignment
        void operator=(const engineTime&);

        //- adjust read time values
        void timeAdjustment();

public:

    // Constructors

        //- Construct from objectRegistry arguments
        engineTime
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const fileName& systemName = "system",
            const fileName& constantName = "constant",
            const fileName& dictName = "engineGeometry"
        );

    //- Destructor
    virtual ~engineTime()
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

            //- Return the engine geometry dictionary
            const dictionary& engineDict() const
            {
                return dict_;
            }

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
            scalar theta() const;

            //- Return current crank-angle translated to a single revolution
            //  (value between -180 and 180 with 0 = top dead centre)
            scalar thetaRevolution() const;

            //- Return crank-angle increment
            scalar deltaTheta() const;

            //- Return current piston position
            dimensionedScalar pistonPosition() const;

            //- Return piston displacement for current time step
            dimensionedScalar pistonDisplacement() const;

            //- Return piston speed for current time step
            dimensionedScalar pistonSpeed() const;


        // Member functions overriding the virtual functions in time

            //- Convert the user-time (CA deg) to real-time (s).
            virtual scalar userTimeToTime(const scalar theta) const;

            //- Convert the real-time (s) into user-time (CA deg)
            virtual scalar timeToUserTime(const scalar t) const;

            //- Read the control dictionary and set the write controls etc.
            virtual void readDict();


        // Edit

            //- Read the controlDict and set all the parameters
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
