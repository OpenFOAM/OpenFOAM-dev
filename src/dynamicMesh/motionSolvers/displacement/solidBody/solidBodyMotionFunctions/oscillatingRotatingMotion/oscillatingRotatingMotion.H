/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::oscillatingRotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    oscillatingRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingRotatingMotion_H
#define oscillatingRotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                  Class oscillatingRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingRotatingMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        vector amplitude_;

        //- Radial velocity
        scalar omega_;


public:

    //- Runtime type information
    TypeName("oscillatingRotatingMotion");


    // Constructors

        //- Construct from components
        oscillatingRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        oscillatingRotatingMotion(const oscillatingRotatingMotion&);

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new oscillatingRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~oscillatingRotatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const oscillatingRotatingMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
