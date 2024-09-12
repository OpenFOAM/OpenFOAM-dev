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
    Foam::Function1s::Shift

Description
    Function1 which shifts a given 'value' function along the x-axis

Usage
    Example:
    \verbatim
        <name>
        {
            type    shift;

            start   10 [CAD];

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
        Symbol | Description        | Data type       | Default
        start  | Start x coordinate | scalar          |
        value  | Function to shift  | Function1<Type> |
    \endvartable

SourceFiles
    Shift.C

\*---------------------------------------------------------------------------*/

#ifndef Shift_H
#define Shift_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Shift Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Shift
:
    public FieldFunction1<Type, Shift<Type>>
{
    // Private Data

        //- Argument offset
        const scalar start_;

        //- Value function
        const autoPtr<Function1<Type>> value_;


public:

    // Runtime type information
    TypeName("shift");


    // Constructors

        //- Construct from name and dictionary
        Shift
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Shift(const Shift<Type>& se);


    //- Destructor
    virtual ~Shift();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Shift<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ShiftI.H"

#ifdef NoRepository
    #include "Shift.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
