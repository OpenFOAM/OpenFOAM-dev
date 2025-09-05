/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::interpolations::cellPoint

Description
    Piecewise-linear interpolation method. Uses volPointInterpolation to create
    values on the points. Then uses a face-diagonal tetrahedral decomposition
    to linearly interpolate between the values at the points and at the
    cell-centres.

SourceFiles
    cellPoint_interpolationI.H
    cellPoint_interpolation.C

\*---------------------------------------------------------------------------*/

#ifndef cellPoint_interpolation_H
#define cellPoint_interpolation_H

#include "interpolation.H"
#include "barycentricTensor.H"
#include "volPointInterpolation_interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolations
{

// Forward declaration of classes
template<class Type>
class cellPoint;

/*---------------------------------------------------------------------------*\
                        Class cellPointBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellPointBase
:
    public fieldInterpolation<Type, cellPoint<Type>>,
    public volPointInterpolation<Type>
{
protected:

    // Protected Member Functions

        //- Return the transform associated with the given tetrahedron
        inline void tetTransform
        (
            const tetIndices& tetIs,
            vector& x0,
            scalar& detA,
            barycentricTensor& T
        ) const;

        //- Return the transform associated with the face-triangle of the given
        //  tetrahedron
        inline void tetFaceTriTransform
        (
            const tetIndices& tetIs,
            vector& x0,
            scalar& detA,
            barycentricTensor& T
        ) const;

        //- Return the barycentric origin for the above transforms
        inline barycentric y0() const;

        //- Find the tetrahedron (or triangle) containing the given point, and
        //  also generate and return the associated transformation
        inline tetIndices findTet
        (
            const point& x,
            const label celli,
            const label facei,
            vector& x0,
            scalar& detA,
            barycentricTensor& T
        ) const;


public:

    // Constructors

        //- Construct from components
        cellPointBase(const VolField<Type>& psi);

        //- Copy construct
        cellPointBase(const cellPointBase<Type>& i);

        //- Construct from components
        cellPointBase(const VolField<Type>& psi, tmp<PointField<Type>> psip);


    // Member Functions

        //- Inherit interpolate from interpolation
        using interpolation<Type>::interpolate;

        //- Interpolate to the given point in the given cell
        inline Type interpolate
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;

        //- Interpolate to the coordinates in the given tetrahedron
        inline Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const;
};


/*---------------------------------------------------------------------------*\
                        Class cellPointGradBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellPointGradBase
:
    public cellPointBase<Type>
{
public:

    // Constructors

        //- Inherit base class constructors
        using cellPointBase<Type>::cellPointBase;


    // Member Functions

        //- Inherit interpolate from interpolation
        using interpolation<Type>::interpolateGrad;

        //- Interpolate the gradient to the given point in the given cell
        inline interpolationGradType<Type> interpolateGrad
        (
            const tetIndices& tetIs,
            const scalar detA,
            const barycentricTensor T
        ) const;

        //- Interpolate the gradient to the given point in the given cell
        inline interpolationGradType<Type> interpolateGrad
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;

        //- Interpolate the gradient to the coordinates in the given tetrahedron
        inline interpolationGradType<Type> interpolateGrad
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const;
};


/*---------------------------------------------------------------------------*\
                 Typedef cellPointBaseOrGradBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
using cellPointBaseOrGradBase =
    typename std::conditional
    <
        (pTraits<Type>::rank > 1),
        cellPointBase<Type>,
        cellPointGradBase<Type>
    >::type;


/*---------------------------------------------------------------------------*\
                          Class cellPoint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellPoint
:
    public cellPointBaseOrGradBase<Type>
{
    // Private Type Definitions

        //- The base class
        typedef cellPointBaseOrGradBase<Type> Base;


public:

    //- Runtime type information
    TypeName("cellPoint");


    // Constructors

        //- Inherit base class constructors
        using Base::Base;

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new cellPoint<Type>(*this)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellPoint_interpolationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellPoint_interpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
