/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::Function2s::OneConstant

Description
    Templated function of two variables that returns the corresponding 1 (one).

    Usage:
    \verbatim
        <name> one;
    \endverbatim

SourceFiles
    OneConstant.C

\*---------------------------------------------------------------------------*/

#ifndef OneConstant2_H
#define OneConstant2_H

#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                         Class OneConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class OneConstant
:
    public FieldFunction2<Type, OneConstant<Type>>
{

public:

    // Runtime type information
    TypeName("one");


    // Constructors

        //- Construct from name
        OneConstant(const word& name);

        //- Construct from name and dictionary
        OneConstant(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~OneConstant();


    // Member Functions

        //- Return One
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Return One field
        virtual inline tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;

        //- Write in dictionary format
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OneConstant<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OneConstant2I.H"

#ifdef NoRepository
    #include "OneConstant2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
