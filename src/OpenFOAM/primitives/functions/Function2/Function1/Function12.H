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
    Foam::Function2s::Function12

Description
    Function2 which returns a Function1 using just one of the arguments given
    to the function2. The function1 is specified as value1 or value2, depending
    on which argument it is to be evaluated with.

    Example for a scalar:
    \verbatim
        <name>
        {
            type    function1;
            value1  constant 10; // <-- The value1 function is evaluated with
                                 //     the first argument
        }
    \endverbatim

    Example for a vector:
    \verbatim
        <name>
        {
            type    function1;
            value2  table       // <-- The value2 function is evaluated with
            (                   //     the second argument
                (0.00 (0 0 0))
                (0.35 (0 0 1))
                (0.71 (0 0 0))
            );
        }
    \endverbatim

SourceFiles
    Function12.C

\*---------------------------------------------------------------------------*/

#ifndef Function12_H
#define Function12_H

#include "Function1.H"
#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class Function12 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function12
:
    public FieldFunction2<Type, Function12<Type>>
{
    // Private Data

        //- The index of the argument to use. 0 or 1, hence a bool.
        bool index_;

        //- Function
        autoPtr<Function1<Type>> f_;


public:

    // Runtime type information
    TypeName("function1");


    // Constructors

        //- Construct from name and dictionary
        Function12
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Function12(const Function12<Type>& se);


    //- Destructor
    virtual ~Function12();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Function12<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Function12I.H"

#ifdef NoRepository
    #include "Function12.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
