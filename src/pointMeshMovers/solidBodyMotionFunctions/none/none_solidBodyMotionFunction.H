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

Class
    Foam::solidBodyMotionFunctions::none

Description
    Dummy motion function for stationary regions

SourceFiles
    none_solidBodyMotionFunction.C

\*---------------------------------------------------------------------------*/

#ifndef none_solidBodyMotionFunction_H
#define none_solidBodyMotionFunction_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                        Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public solidBodyMotionFunction
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        none
        (
            const word& name,
            const Time& runTime
        );

        //- Construct from components
        none
        (
            const word& name,
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new none(name_, time_)
            );
        }


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Return the static solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Dummy read
        virtual bool read(const dictionary&);

        //- Dummy write
        virtual void writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
