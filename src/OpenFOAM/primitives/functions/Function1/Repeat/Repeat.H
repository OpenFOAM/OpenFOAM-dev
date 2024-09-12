/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::Function1s::Repeat

Description
    Function1 which repeats a given 'value' function with a given period or
    frequency, and with an optional shift along the x-axis

Usage
    Example:
    \verbatim
        <name>
        {
            type    repeat;

            period  360 [CAD];

            //frequency 0.00277778 [CAD^-1]; // <-- alternative to specifying
                                             //     the period

            start   10 [CAD];                // <-- optional shift

            value
            {
                type    table;
                file    "constant/timeVsFlowRate.foam";
                units   ([ms] [m^3/s]);
            }
        }
    \endverbatim

    Where:
    \vartable
        Symbol    | Description             | Data type       | Default
        period    | Period of repetition    | scalar          |
        frequency | Frequency of repetition | scalar          |
        start     | Start x coordinate      | scalar          | 0
        value     | Function to repeat      | Function1<Type> |
    \endvartable

SourceFiles
    Repeat.C

\*---------------------------------------------------------------------------*/

#ifndef Repeat_H
#define Repeat_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Repeat Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Repeat
:
    public FieldFunction1<Type, Repeat<Type>>
{
    // Private Data

        //- Period of repetition
        const scalar period_;

        //- Argument offset
        const scalar start_;

        //- Value function
        const autoPtr<Function1<Type>> value_;


    // Private Member Functions

        //- Read the period
        static scalar readPeriod
        (
            const unitConversions& units,
            const dictionary& dict
        );


public:

    // Runtime type information
    TypeName("repeat");


    // Constructors

        //- Construct from name and dictionary
        Repeat
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Repeat(const Repeat<Type>& se);


    //- Destructor
    virtual ~Repeat();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Repeat<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RepeatI.H"

#ifdef NoRepository
    #include "Repeat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
