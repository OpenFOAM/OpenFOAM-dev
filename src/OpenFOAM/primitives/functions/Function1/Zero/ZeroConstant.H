/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::Function1Types::ZeroConstant

Description
    Templated function that returns the corresponding 0 (zero).

    Usage:
    \verbatim
        <entryName> zero;
    \endverbatim

SourceFiles
    ZeroConstant.C

\*---------------------------------------------------------------------------*/

#ifndef ZeroConstant_H
#define ZeroConstant_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class ZeroConstant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ZeroConstant
:
    public Function1<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ZeroConstant<Type>&);


public:

    // Runtime type information
    TypeName("zero");


    // Constructors

        //- Construct from entry name and dictionary
        ZeroConstant(const word& entryName, const dictionary& dict);


    //- Destructor
    virtual ~ZeroConstant();


    // Member Functions

        //- Return constant value
        virtual inline Type value(const scalar) const;

        //- Integrate between two values
        virtual inline Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ZeroConstantI.H"

#ifdef NoRepository
    #include "ZeroConstant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
