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
    Foam::Function2s::Radial

Description
    Function2 which returns a Function1 of the magnitude of the two-dimensional
    vector with components equal to the input arguments.

    Example
    \verbatim
        <name>
        {
            type    radial;
            value   table
            (
                (0.00 (0 0 0))
                (0.35 (0 0 1))
                (0.71 (0 0 0))
            );
        }
    \endverbatim

SourceFiles
    Radial.C

\*---------------------------------------------------------------------------*/

#ifndef Radial2_H
#define Radial2_H

#include "Function1.H"
#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                            Class Radial Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Radial
:
    public FieldFunction2<Type, Radial<Type>>
{
    // Private Data

        //- Function of the radius
        autoPtr<Function1<Type>> value_;


public:

    // Runtime type information
    TypeName("radial");


    // Constructors

        //- Construct from name and dictionary
        Radial
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Radial(const Radial<Type>& se);


    //- Destructor
    virtual ~Radial();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Radial<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Radial2I.H"

#ifdef NoRepository
    #include "Radial2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
