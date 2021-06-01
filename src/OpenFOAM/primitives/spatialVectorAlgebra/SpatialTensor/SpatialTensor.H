/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::SpatialTensor

Description
    Templated 3D spatial tensor derived from MatrixSpace used to represent
    transformations of spatial vectors and the angular and linear inertia of
    rigid bodies.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    SpatialTensorI.H

See also
    Foam::MatrixSpace
    Foam::Tensor
    Foam::SpatialVector

\*---------------------------------------------------------------------------*/

#ifndef SpatialTensor_H
#define SpatialTensor_H

#include "Tensor.H"
#include "SpatialVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SpatialTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SpatialTensor
:
    public MatrixSpace<SpatialTensor<Cmpt>, Cmpt, 6, 6>
{

public:

    // Member constants

        //- Rank of Tensor is 2
        static const direction rank = 2;


    // Static Data Members

        //- Identity matrix for square matrices
        static const SpatialTensor I;


    // Constructors

        //- Construct null
        inline SpatialTensor();

        //- Construct initialised to zero
        inline SpatialTensor(const Foam::zero);

        //- Construct given MatrixSpace of the same rank
        inline SpatialTensor(const typename SpatialTensor::msType&);

        //- Construct given 4 tensor blocks
        inline SpatialTensor
        (
            const Tensor<Cmpt>& t00, const Tensor<Cmpt>& t01,
            const Tensor<Cmpt>& t10, const Tensor<Cmpt>& t11
        );

        //- Construct given 36 components
        inline SpatialTensor
        (
            const Cmpt& t00, const Cmpt& t01, const Cmpt& t02,
            const Cmpt& t03, const Cmpt& t04, const Cmpt& t05,

            const Cmpt& t10, const Cmpt& t11, const Cmpt& t12,
            const Cmpt& t13, const Cmpt& t14, const Cmpt& t15,

            const Cmpt& t20, const Cmpt& t21, const Cmpt& t22,
            const Cmpt& t23, const Cmpt& t24, const Cmpt& t25,

            const Cmpt& t30, const Cmpt& t31, const Cmpt& t32,
            const Cmpt& t33, const Cmpt& t34, const Cmpt& t35,

            const Cmpt& t40, const Cmpt& t41, const Cmpt& t42,
            const Cmpt& t43, const Cmpt& t44, const Cmpt& t45,

            const Cmpt& t50, const Cmpt& t51, const Cmpt& t52,
            const Cmpt& t53, const Cmpt& t54, const Cmpt& t55
        );

        //- Construct from Istream
        inline SpatialTensor(Istream&);
};


template<class Cmpt>
class typeOfTranspose<Cmpt, SpatialTensor<Cmpt>>
{
public:

    typedef SpatialTensor<Cmpt> type;
};


template<class Cmpt>
class typeOfOuterProduct<Cmpt, SpatialVector<Cmpt>, SpatialVector<Cmpt>>
{
public:

    typedef SpatialTensor<Cmpt> type;
};


template<class Cmpt>
class typeOfInnerProduct<Cmpt, SpatialTensor<Cmpt>, SpatialVector<Cmpt>>
{
public:

    typedef SpatialVector<Cmpt> type;
};


template<class Cmpt>
class typeOfInnerProduct<Cmpt, SpatialTensor<Cmpt>, SpatialTensor<Cmpt>>
{
public:

    typedef SpatialTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SpatialTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
