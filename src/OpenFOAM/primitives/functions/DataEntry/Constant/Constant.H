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
    Foam::Constant

Description
    Templated basic entry that holds a constant value.

    Usage - for entry \<entryName\> having the value <value>:
    \verbatim
        <entryName>    constant  <value>
    \endverbatim

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "DataEntry.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Constant;

template<class Type>
Ostream& operator<<(Ostream&, const Constant<Type>&);

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public DataEntry<Type>
{
    // Private data

        //- Constant value
        Type value_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Constant<Type>&);


public:

    // Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from entry name and Istream
        Constant(const word& entryName, const dictionary& dict);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);

        //- Construct and return a clone
        virtual tmp<DataEntry<Type> > clone() const
        {
            return tmp<DataEntry<Type> >(new Constant<Type>(*this));
        }


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant value
        Type value(const scalar) const;

        //- Integrate between two values
        Type integrate(const scalar x1, const scalar x2) const;

        //- Return dimensioned constant value
        dimensioned<Type> dimValue(const scalar) const;

        //- Integrate between two values and return dimensioned type
        dimensioned<Type> dimIntegrate
        (
            const scalar x1,
            const scalar x2
        ) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Constant<Type>& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Constant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
