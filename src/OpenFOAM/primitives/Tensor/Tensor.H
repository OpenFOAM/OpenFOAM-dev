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
    Foam::Tensor

Description
    Templated 3D tensor derived from MatrixSpace adding construction from
    9 components, element access using xx(), xy() etc. member functions and
    the inner-product (dot-product) and outer-product of two Vectors
    (tensor-product) operators.

SourceFiles
    TensorI.H

See also
    Foam::MatrixSpace
    Foam::Vector

\*---------------------------------------------------------------------------*/

#ifndef Tensor_H
#define Tensor_H

#include "MatrixSpace.H"
#include "Vector.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Cmpt>
class SymmTensor;

template<class Cmpt>
class DiagTensor;

/*---------------------------------------------------------------------------*\
                           Class Tensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor
:
    public MatrixSpace<Tensor<Cmpt>, Cmpt, 3, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Tensor<label> labelType;


    // Member constants

        //- Rank of Tensor is 2
        static const direction rank = 2;


    // Static Data Members

        static const Tensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };


    // Constructors

        //- Construct null
        inline Tensor();

        //- Construct initialised to zero
        inline Tensor(const Foam::zero);

        //- Construct given MatrixSpace of the same rank
        template<class Cmpt2>
        inline Tensor(const MatrixSpace<Tensor<Cmpt2>, Cmpt2, 3, 3>&);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline Tensor(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>&);

        //- Construct given SphericalTensor
        inline Tensor(const SphericalTensor<Cmpt>&);

        //- Construct given SymmTensor
        inline Tensor(const SymmTensor<Cmpt>&);

        //- Construct given DiagTensor
        inline Tensor(const DiagTensor<Cmpt>&);

        //- Construct given triad
        inline Tensor(const Vector<Vector<Cmpt>>&);

        //- Construct given the three vector components
        inline Tensor
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z
        );

        //- Construct given the nine components
        inline Tensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
            const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
            const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
        );

        //- Construct from a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        Tensor
        (
            const Block2<Tensor<Cmpt>, BRowStart, BColStart>& block
        );

        //- Construct from Istream
        inline Tensor(Istream&);


    // Member Functions

        // Component access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& xz() const;
            inline const Cmpt& yx() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& yz() const;
            inline const Cmpt& zx() const;
            inline const Cmpt& zy() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& xz();
            inline Cmpt& yx();
            inline Cmpt& yy();
            inline Cmpt& yz();
            inline Cmpt& zx();
            inline Cmpt& zy();
            inline Cmpt& zz();

        // Row-vector access.

            inline Vector<Cmpt> x() const;
            inline Vector<Cmpt> y() const;
            inline Vector<Cmpt> z() const;
            inline Vector<Cmpt> vectorComponent(const direction) const;

        //- Return transpose
        inline Tensor<Cmpt> T() const;

        //- Return inverse
        inline Tensor<Cmpt> inv() const;


    // Member Operators

        //- Inner-product with a Tensor
        inline void operator&=(const Tensor<Cmpt>&);

        //- Inherit MatrixSpace assignment operators
        using Tensor::msType::operator=;

        //- Assign to an equivalent vector space
        template<class Cmpt2>
        inline void operator=(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>&);

        //- Assign to a SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);

        //- Assign to a SymmTensor
        inline void operator=(const SymmTensor<Cmpt>&);

        //- Assign to a DiagTensor
        inline void operator=(const DiagTensor<Cmpt>&);

        //- Assign to a triad
        inline void operator=(const Vector<Vector<Cmpt>>&);
};


template<class Cmpt>
class typeOfRank<Cmpt, 2>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfTranspose<Cmpt, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "TensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
