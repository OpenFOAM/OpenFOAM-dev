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
    Foam::Function2s::Product

Description
    Function2 which returns the product of two independent Function1-s of the
    two input arguments. The two function1s are specified as value1 and value2.

    Example to scale a table of vectors in the first argument with a ramp in
    the second argument:
    \verbatim
        <name>
        {
            type            product;
            value1<vector>  table
            (
                (0.00 (0 0 0))
                (0.25 (1 0 0))
                (0.50 (0 0 0))
            );
            value2<scalar>
            {
                type        linearRamp;
                start       1;
                duration    4;
            }
        }
    \endverbatim

    Note that only one type specification (the \<vector\>/\<scalar\>/... part)
    is needed in general for the value entries, and no type specifications
    are needed if the function is scalar.

SourceFiles
    Product.C

\*---------------------------------------------------------------------------*/

#ifndef Product2_H
#define Product2_H

#include "Function1.H"
#include "Function2.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{

/*---------------------------------------------------------------------------*\
                  Class ProductValueTypeIsValid Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class ValueType, class Enable = void>
struct ProductValueTypeIsValid
{
    static const bool value = false;
};


template<class Type>
struct ProductValueTypeIsValid<Type, scalar>
{
    static const bool value = true;
};


template<class Type, class ValueType>
struct ProductValueTypeIsValid
<
    Type,
    ValueType,
    typename std::enable_if
    <
        (pTraits<ValueType>::rank > 0)
     && (pTraits<ValueType>::rank <= pTraits<Type>::rank)
    >::type
>
{
    static const bool value =
        std::is_same
        <
            Type,
            typename outerProduct
            <
                ValueType,
                typename typeOfRank
                <
                    scalar,
                    pTraits<Type>::rank - pTraits<ValueType>::rank
                >::type
            >::type
        >::value;
};


/*---------------------------------------------------------------------------*\
                       Class ProductFilter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct ProductFilter
{
    inline const Type& operator()(const Type& value) const
    {
        return value;
    }

    template<class WrongType>
    const Type& operator()(const WrongType& value) const
    {
        return NullObjectRef<Type>();
    }
};


/*---------------------------------------------------------------------------*\
                      Class ProductFunction1s Declaration
\*---------------------------------------------------------------------------*/

template<class Type, direction rank = pTraits<Type>::rank>
class ProductFunction1s
:
    public ProductFunction1s<Type, rank - 1>
{
    // Private Typedefs

        //- Type of the function
        typedef
            Function1
            <
                typename std::conditional
                <
                    rank == pTraits<Type>::rank,
                    Type,
                    typename typeOfRank<scalar, rank>::type
                >::type
            > function1Type;


public:

    // Public Data

        //- Functions
        Pair<autoPtr<function1Type>> fs;


    // Constructors

        //- Construct from a dictionary
        ProductFunction1s
        (
            const unitConversions& units,
            const dictionary& dict,
            const Pair<Tuple2<word, label>>& typeAndRanks
        );

        //- Copy construct
        ProductFunction1s(const ProductFunction1s<Type, rank>& p2f1s);


    // Public Member Functions

        //- Evaluate
        inline Type value(const scalar x, const scalar y) const;

        //- Write to a stream
        void write(Ostream& os, const unitConversions& units) const;
};


/*---------------------------------------------------------------------------*\
                    Class ProductFunction1s<0> Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ProductFunction1s<Type, 0>
{
public:

    // Public Data

        //- Functions
        Pair<autoPtr<Function1<scalar>>> fs;


    // Constructors

        //- Construct from a dictionary
        ProductFunction1s
        (
            const unitConversions& units,
            const dictionary& dict,
            const Pair<Tuple2<word, label>>& typeAndRanks
        );

        //- Copy construct
        ProductFunction1s(const ProductFunction1s<Type, 0>& p2f1s);


    // Public Member Functions

        //- Evaluate
        inline Type value(const scalar x, const scalar y) const;

        //- Write to a stream
        void write(Ostream& os, const unitConversions& units) const;
};


/*---------------------------------------------------------------------------*\
                            Class Product Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Product
:
    public FieldFunction2<Type, Product<Type>>
{
    // Private Data

        //- Functions
        ProductFunction1s<Type> fs_;


public:

    // Runtime type information
    TypeName("product");


    // Constructors

        //- Construct from name and dictionary
        Product
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        Product(const Product<Type>& se);


    //- Destructor
    virtual ~Product();


    // Member Functions

        //- Return value
        virtual inline Type value(const scalar x, const scalar y) const;

        //- Write in dictionary format
        virtual void write(Ostream& os, const unitConversions& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Product<Type>&) = delete;
};


// Global Functions

//- Return the name of the value entry for the given argument
inline word valueName(const direction argument);

//- Return the name of the value entry for the given argument and type name
inline word valueName(const direction argument, const word& typeName);

//- Return the name of the value entry for the given argument and type
template<class Type>
inline word valueName(const direction argument);

//- Return the name of the value entry for the given argument and type and rank
inline word valueName
(
    const direction argument,
    const Tuple2<word, label>& typeAndRank
);

//- Lookup the type and rank for the value entry for the given argument
template<class Type, class ValueType>
void lookupValueTypeAndRank
(
    const dictionary& dict,
    const direction argument,
    Tuple2<word, label>& typeAndRank,
    label& found
);

//- Lookup the type and rank for the value entry for the given argument
template<class Type>
Tuple2<word, label> lookupValueTypeAndRank
(
    const dictionary& dict,
    const direction argument
);

//- Lookup the types and ranks for the value entries
template<class Type>
Pair<Tuple2<word, label>> lookupValueTypeAndRanks(const dictionary& dict);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function2s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Product2I.H"

#ifdef NoRepository
    #include "Product2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
