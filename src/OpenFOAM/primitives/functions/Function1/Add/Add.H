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
    Foam::Function1s::Add

Description
    Function1 which adds the value of two Function1-s.

Usage
    Example, to add an offset to a table:
    \verbatim
        <name>
        {
            type      add;

            value1
            {
                type            table;
                format          csv;
                nHeaderLine     1;
                units           ([CAD] [kPa]);
                columns         (0 1);
                mergeSeparators no;
                file            "constant/expData/pInlet";
            }

            value2      constant 1 [bar];
        }
    \endverbatim

    Where:
    \table
        Property | Description                                    | Required
        value1   | Function of type Function1<Type>               | yes
        value2   | Function of type Function1<Type>               | yes
    \endtable

SourceFiles
    Add.C

\*---------------------------------------------------------------------------*/

#ifndef Add_H
#define Add_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                            Class Add Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Add
:
    public FieldFunction1<Type, Add<Type>>
{
    // Private Data

        //- Value1 function
        const autoPtr<Function1<Type>> value1_;

        //- Value2 function
        const autoPtr<Function1<Type>> value2_;


public:

    // Runtime type information
    TypeName("add");


    // Constructors

        //- Construct from name and dictionary
        Add
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Add(const Add<Type>& se);


    //- Destructor
    virtual ~Add();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Is this function guaranteed to be constant?
        virtual inline bool constant() const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Add<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AddI.H"

#ifdef NoRepository
    #include "Add.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
