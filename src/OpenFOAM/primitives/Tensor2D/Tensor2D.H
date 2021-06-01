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
    Foam::Tensor2D

Description
    Templated 2D tensor derived from VectorSpace adding construction from
    4 components, element access using xx(), xy(), yx() and yy() member
    functions and the iner-product (dot-product) and outer-product of two
    Vector2Ds (tensor-product) operators.

SourceFiles
    Tensor2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Tensor2D_H
#define Tensor2D_H

#include "Vector2D.H"
#include "SphericalTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Cmpt>
class SymmTensor2D;

/*---------------------------------------------------------------------------*\
                          Class Tensor2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor2D
:
    public VectorSpace<Tensor2D<Cmpt>, Cmpt, 4>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Tensor2D<label> labelType;


    // Member constants

        //- Rank of Tensor2D is 2
        static const direction rank = 2;


    // Static Data Members

        static const Tensor2D I;


    //- Component labeling enumeration
    enum components { XX, XY, YX, YY };


    // Constructors

        //- Construct null
        inline Tensor2D();

        //- Construct initialised to zero
        inline Tensor2D(const Foam::zero);

        //- Construct given VectorSpace
        inline Tensor2D(const VectorSpace<Tensor2D<Cmpt>, Cmpt, 4>&);

        //- Construct given SymmTensor2D
        inline Tensor2D(const SymmTensor2D<Cmpt>&);

        //- Construct given SphericalTensor2D
        inline Tensor2D(const SphericalTensor2D<Cmpt>&);

        //- Construct given the two vectors
        inline Tensor2D
        (
            const Vector2D<Cmpt>& x,
            const Vector2D<Cmpt>& y
        );

        //- Construct given the four components
        inline Tensor2D
        (
            const Cmpt txx, const Cmpt txy,
            const Cmpt tyx, const Cmpt tyy
        );

        //- Construct from Istream
        Tensor2D(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& yx() const;
            inline const Cmpt& yy() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& yx();
            inline Cmpt& yy();

            // Access vector components.

            inline Vector2D<Cmpt> x() const;
            inline Vector2D<Cmpt> y() const;

        //- Transpose
        inline Tensor2D<Cmpt> T() const;


    // Member Operators

        //- Copy SymmTensor2D
        inline void operator=(const SymmTensor2D<Cmpt>&);

        //- Copy SphericalTensor2D
        inline void operator=(const SphericalTensor2D<Cmpt>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "Tensor2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
