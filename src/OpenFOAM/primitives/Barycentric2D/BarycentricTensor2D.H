/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::BarycentricTensor2D

Description
    Templated 3x3 tensor derived from VectorSpace. Has 9 components. Can
    represent a barycentric transformation as a matrix-barycentric inner-
    product. Can alternatively represent an inverse barycentric transformation
    as a vector-matrix inner-product.

SourceFiles
    BarycentricTensor2DI.H

\*---------------------------------------------------------------------------*/

#ifndef BarycentricTensor2D_H
#define BarycentricTensor2D_H

#include "Barycentric2D.H"
#include "Tensor2D.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BarycentricTensor2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class BarycentricTensor2D
:
    public MatrixSpace<BarycentricTensor2D<Cmpt>, Cmpt, 3, 3>
{
public:

    //- Equivalent type of labels used for valid component indexing
    typedef Tensor2D<label> labelType;


    // Member constants

        //- Rank of BarycentricTensor2D is 2
        static const direction rank = 2;


    //- Component labeling enumeration
    enum components { XA, XB, XC, YA, YB, YC, ZA, ZB, ZC };


    // Constructors

        //- Construct null
        BarycentricTensor2D();

        //- Construct initialised to zero
        BarycentricTensor2D(const Foam::zero);

        //- Construct given three barycentric components (rows)
        BarycentricTensor2D
        (
            const Barycentric2D<Cmpt>& x,
            const Barycentric2D<Cmpt>& y,
            const Barycentric2D<Cmpt>& z
        );

        //- Construct given three vector components (columns)
        BarycentricTensor2D
        (
            const Vector<Cmpt>& a,
            const Vector<Cmpt>& b,
            const Vector<Cmpt>& c
        );


    // Member Functions

        // Row-barycentric access

            inline Barycentric2D<Cmpt> x() const;
            inline Barycentric2D<Cmpt> y() const;
            inline Barycentric2D<Cmpt> z() const;

        // Column-vector access

            inline Vector<Cmpt> a() const;
            inline Vector<Cmpt> b() const;
            inline Vector<Cmpt> c() const;
};


template<class Cmpt>
class typeOfTranspose<Cmpt, BarycentricTensor2D<Cmpt>>
{
public:

    typedef void type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BarycentricTensor2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
