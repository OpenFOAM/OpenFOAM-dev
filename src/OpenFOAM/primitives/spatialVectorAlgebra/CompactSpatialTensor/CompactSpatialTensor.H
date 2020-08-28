/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::CompactSpatialTensor

Description
    Templated 3D compact spatial tensor derived from MatrixSpace used to
    represent transformations of spatial vectors and the angular and linear
    inertia of rigid bodies.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    CompactSpatialTensorI.H

See also
    Foam::MatrixSpace
    Foam::SpatialTensor

\*---------------------------------------------------------------------------*/

#ifndef CompactSpatialTensor_H
#define CompactSpatialTensor_H

#include "SpatialTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class CompactSpatialTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class CompactSpatialTensor
:
    public MatrixSpace<CompactSpatialTensor<Cmpt>, Cmpt, 6, 3>
{

public:

    // Constructors

        //- Construct null
        inline CompactSpatialTensor();

        inline CompactSpatialTensor(const Foam::zero);

        //- Construct given MatrixSpace of the same rank
        inline CompactSpatialTensor
        (
            const typename CompactSpatialTensor::msType&
        );

        //- Construct given 18 components
        inline CompactSpatialTensor
        (
            const Cmpt& t00, const Cmpt& t01, const Cmpt& t02,
            const Cmpt& t10, const Cmpt& t11, const Cmpt& t12,
            const Cmpt& t20, const Cmpt& t21, const Cmpt& t22,
            const Cmpt& t30, const Cmpt& t31, const Cmpt& t32,
            const Cmpt& t40, const Cmpt& t41, const Cmpt& t42,
            const Cmpt& t50, const Cmpt& t51, const Cmpt& t52
        );

        //- Construct from Istream
        inline CompactSpatialTensor(Istream&);
};


template<class Cmpt>
class typeOfInnerProduct<Cmpt, CompactSpatialTensor<Cmpt>, Tensor<Cmpt>>
{
public:

    typedef CompactSpatialTensor<Cmpt> type;
};


template<class Cmpt>
class typeOfInnerProduct<Cmpt, CompactSpatialTensor<Cmpt>, Vector<Cmpt>>
{
public:

    typedef SpatialVector<Cmpt> type;
};


template<class Cmpt>
class typeOfInnerProduct<Cmpt, SpatialTensor<Cmpt>, CompactSpatialTensor<Cmpt>>
{
public:

    typedef CompactSpatialTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "CompactSpatialTensorT.H"
#include "CompactSpatialTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
