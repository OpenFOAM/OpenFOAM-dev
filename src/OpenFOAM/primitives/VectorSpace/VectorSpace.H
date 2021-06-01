/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::VectorSpace

Description
    Templated vector space.

    Template arguments are the Form the vector space will be used to create,
    the type of the elements and the number of elements.

SourceFiles
    VectorSpaceI.H
    VectorSpace.C

\*---------------------------------------------------------------------------*/

#ifndef VectorSpace_H
#define VectorSpace_H

#include "direction.H"
#include "scalar.H"
#include "word.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Form, class Cmpt, direction Ncmpts> class VectorSpace;

template<class Form, class Cmpt, direction Ncmpts>
void writeEntry(Ostream& os, const VectorSpace<Form, Cmpt, Ncmpts>& value);

template<class Form, class Cmpt, direction Ncmpts>
Istream& operator>>
(
    Istream&,
    VectorSpace<Form, Cmpt, Ncmpts>&
);

template<class Form, class Cmpt, direction Ncmpts>
Ostream& operator<<
(
    Ostream&,
    const VectorSpace<Form, Cmpt, Ncmpts>&
);


/*---------------------------------------------------------------------------*\
                         Class VectorSpace Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Cmpt, direction Ncmpts>
class VectorSpace
{

public:

    //- The components of this vector space
    Cmpt v_[Ncmpts];

    //- VectorSpace type
    typedef VectorSpace<Form, Cmpt, Ncmpts> vsType;

    //- Component type
    typedef Cmpt cmptType;


    // Static constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Number of components in this vector space
        static const direction nComponents = Ncmpts;


        // VectorSpace currently defaults to a column-vector
        // This will be removed when column-vector is introduced
        // as a specialisation
        static const direction mRows = Ncmpts;
        static const direction nCols = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const Form zero;
        static const Form one;
        static const Form max;
        static const Form min;
        static const Form rootMax;
        static const Form rootMin;


    // Sub-Block Classes

        //- Const sub-block type
        template
        <
            class SubVector,
            direction BStart
        >
        class ConstBlock
        {
            const vsType& vs_;

        public:

            //- Number of components in this vector space
            static const direction nComponents = SubVector::nComponents;

            //- Construct for a given vector
            inline ConstBlock(const vsType& vs);

            //- [i] const element access operator
            inline const Cmpt& operator[](const direction i) const;

            //- (i, 0) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction
            ) const;
        };


    // Constructors

        //- Construct null
        inline VectorSpace();

        //- Construct initialised to zero
        inline VectorSpace(const Foam::zero);

        //- Construct from Istream
        VectorSpace(Istream&);

        //- Construct as copy of a VectorSpace with the same size
        template<class Form2, class Cmpt2>
        inline explicit VectorSpace(const VectorSpace<Form2, Cmpt2, Ncmpts>&);


    // Member Functions

        //- Return the number of elements in the VectorSpace = Ncmpts.
        inline static direction size();

        inline const Cmpt& component(const direction) const;
        inline Cmpt& component(const direction);

        inline void component(Cmpt&, const direction) const;
        inline void replace(const direction, const Cmpt&);

        //- Return a VectorSpace with all elements = s
        inline static Form uniform(const Cmpt& s);

        template<class SubVector, direction BStart>
        inline const ConstBlock<SubVector, BStart> block() const;


    // Member Operators

        inline const Cmpt& operator[](const direction) const;
        inline Cmpt& operator[](const direction);

        inline void operator+=(const VectorSpace<Form, Cmpt, Ncmpts>&);
        inline void operator-=(const VectorSpace<Form, Cmpt, Ncmpts>&);

        inline void operator=(const Foam::zero);
        inline void operator*=(const scalar);
        inline void operator/=(const scalar);


    // IOstream Operators

        friend Istream& operator>> <Form, Cmpt, Ncmpts>
        (
            Istream&,
            VectorSpace<Form, Cmpt, Ncmpts>&
        );

        friend Ostream& operator<< <Form, Cmpt, Ncmpts>
        (
            Ostream&,
            const VectorSpace<Form, Cmpt, Ncmpts>&
        );
};


// * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * * //

//- Return a string representation of a VectorSpace
template<class Form, class Cmpt, direction Ncmpts>
word name(const VectorSpace<Form, Cmpt, Ncmpts>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VectorSpaceI.H"

#ifdef NoRepository
    #include "VectorSpace.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
