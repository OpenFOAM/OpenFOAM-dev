/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::Function1s::Polynomial

Description
    Arbitrary order polynomial Function1.

    The coefficients in the polynomial are specified in order of increasing
    exponent, starting from the constant coefficient (i.e., zero exponent).

    For example, the polynomial y = -1 + x^2 + 2x^3 can be specified as:

    \verbatim
        <name>  polynomial (-1 0 1 2);
    \endverbatim

    Or, alternatively:

    \verbatim
        <name>
        {
            type    polynomial;
            coeffs  (-1 0 1 2);
        }
    \endverbatim

SourceFiles
    Polynomial1.C

\*---------------------------------------------------------------------------*/

#ifndef Polynomial1_H
#define Polynomial1_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class Polynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Polynomial
:
    public FieldFunction1<Type, Polynomial<Type>>
{
    // Private Data

        //- Polynomial coefficients
        const List<Type> coeffs_;


    // Private Member Functions

        //- Do unit conversions on input data
        static List<Type> convertRead(const unitConversions&, List<Type>);

        //- Do unit conversions in order to output data
        static List<Type> convertWrite(const unitConversions&, List<Type>);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from name and dictionary
        Polynomial
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Construct from name and Istream
        Polynomial
        (
            const word& name,
            const unitConversions& units,
            Istream& is
        );

        //- Copy constructor
        Polynomial(const Polynomial& poly);


    //- Destructor
    virtual ~Polynomial();


    // Member Functions

        // Evaluation

            //- Return Polynomial value as a function of scalar x
            virtual Type value(const scalar x) const;

            //- Integrate between two scalar fields
            virtual Type integral(const scalar x1, const scalar x2) const;


        //- Write in dictionary format
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Polynomial<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Polynomial1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
