/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::PolynomialEntry

Description
    PolynomialEntry container data entry for scalars. Items are stored in a
    list of Tuple2's. Data is input in the form,
    e.g. for an entry \<entryName\> that describes y = x^2 + 2x^3

    \verbatim
        <entryName>  polynomial  [0 0 1 0 0] // optional dimensions
        (
            (1    2)
            (2    3)
        );
    \endverbatim

SourceFiles
    PolynomialEntry.C

\*---------------------------------------------------------------------------*/

#ifndef PolynomialEntry_H
#define PolynomialEntry_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "dimensionSet.H"
#include "DataEntryFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class PolynomialEntry;

// Forward declaration of friend functions
template<class Type>
Ostream& operator<<(Ostream&, const PolynomialEntry<Type>&);

/*---------------------------------------------------------------------------*\
                        Class PolynomialEntry Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PolynomialEntry
:
    public DataEntry<Type>
{
    // Private data

        //- PolynomialEntry coefficients - list of prefactor, exponent
        List<Tuple2<Type, Type> > coeffs_;

        //- Flag to indicate whether poly can be integrated
        bool canIntegrate_;

        //- The dimension set
        dimensionSet dimensions_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const PolynomialEntry<Type>&);


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        PolynomialEntry(const word& entryName, const dictionary& dict);

        //- Construct from components
        PolynomialEntry
        (
            const word& entryName,
            const List<Tuple2<Type, Type> >&
        );

        //- Copy constructor
        PolynomialEntry(const PolynomialEntry& poly);

        //- Construct and return a clone
        virtual tmp<DataEntry<Type> > clone() const
        {
            return tmp<DataEntry<Type> >(new PolynomialEntry(*this));
        }


    //- Destructor
    virtual ~PolynomialEntry();


    // Member Functions

        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return PolynomialEntry value
            Type value(const scalar x) const;

            //- Integrate between two (scalar) values
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
            const PolynomialEntry<Type>& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PolynomialEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
