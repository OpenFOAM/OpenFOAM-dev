/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::polynomial

Description
    Polynomial container data entry for scalars. Items are stored in a list of
    Tuple2's. Data is input in the form, e.g. for an entry \<entryName\> that
    describes y = x^2 + 2x^3

    \verbatim
        <entryName>   polynomial  [0 0 1 0 0] // optional dimensions
        (
            (1    2)
            (2    3)
        );
    \endverbatim

SourceFiles
    polynomial.C

\*---------------------------------------------------------------------------*/

#ifndef polynomial_H
#define polynomial_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polynomial;

// Forward declaration of friend functions
Ostream& operator<<
(
    Ostream&,
    const polynomial&
);

/*---------------------------------------------------------------------------*\
                        Class polynomial Declaration
\*---------------------------------------------------------------------------*/

class polynomial
:
    public scalarDataEntry
{
    // Private data

        //- Polynomial coefficients - list of prefactor, exponent
        List<Tuple2<scalar, scalar> > coeffs_;

        //- Flag to indicate whether poly can be integrated
        bool canIntegrate_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const polynomial&);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        polynomial(const word& entryName, const dictionary& dict);

        //- Construct from components
        polynomial(const word& entryName, const List<Tuple2<scalar, scalar> >&);

        //- Copy constructor
        polynomial(const polynomial& poly);

        //- Construct and return a clone
        virtual tmp<scalarDataEntry> clone() const
        {
            return tmp<scalarDataEntry>(new polynomial(*this));
        }


    //- Destructor
    virtual ~polynomial();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return polynomial value
            scalar value(const scalar x) const;

            //- Integrate between two (scalar) values
            scalar integrate(const scalar x1, const scalar x2) const;

            //- Return dimensioned constant value
            dimensioned<scalar> dimValue(const scalar) const;

            //- Integrate between two values and return dimensioned type
            dimensioned<scalar> dimIntegrate
            (
                const scalar x1,
                const scalar x2
            ) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const polynomial& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
