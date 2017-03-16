/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::Function1Types::Constant

Description
    Templated function that returns a constant value.

    Usage - for entry \<entryName\> returning the value <value>:
    \verbatim
        <entryName>    constant  <value>
    \endverbatim

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public Function1<Type>
{
    // Private data

        //- Constant value
        Type value_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Constant<Type>&);


public:

    // Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from entry name and value
        Constant(const word& entryName, const Type& val);

        //- Construct from entry name and dictionary
        Constant(const word& entryName, const dictionary& dict);

        //- Construct from entry name and Istream
        //  Reads the constant value without the Function1 type
        //  for backward compatibility
        Constant(const word& entryName, Istream& is);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new Constant<Type>(*this));
        }


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant value
        Type value(const scalar) const;

        //- Integrate between two values
        Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Constant.C"
    #include "Function1New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
