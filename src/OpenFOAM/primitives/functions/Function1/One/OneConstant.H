/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::Function1s::OneConstant

Description
    Templated function that returns the corresponding 1 (one).

    Usage:
    \verbatim
        <name> one;
    \endverbatim

SourceFiles
    OneConstant.C

\*---------------------------------------------------------------------------*/

#ifndef OneConstant_H
#define OneConstant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class OneConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class OneConstant
:
    public FieldFunction1<Type, OneConstant<Type>>
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

        //- Return constant value
        virtual inline Type value(const scalar) const;

        //- Integrate between two values
        virtual inline Type integral(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OneConstant<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OneConstantI.H"

#ifdef NoRepository
    #include "OneConstant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
