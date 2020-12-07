/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

    The coefficients and exponents of the terms in the polynomial are specified
    as a list of Tuple2's e.g. for the polynomial y = x^2 + 2x^3

    \verbatim
        <name>  polynomial
        (
            (1    2)
            (2    3)
        );
    \endverbatim

SourceFiles
    Polynomial1.C

\*---------------------------------------------------------------------------*/

#ifndef Polynomial1_H
#define Polynomial1_H

#include "Function1.H"
#include "Tuple2.H"

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

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<Type, Type>> coeffs_;

        //- Flag to indicate whether poly can be integrald
        bool canIntegrate_;


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from name and dictionary
        Polynomial(const word& name, const dictionary& dict);

        //- Construct from components
        Polynomial
        (
            const word& name,
            const List<Tuple2<Type, Type>>&
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
        virtual void write(Ostream& os) const;


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
