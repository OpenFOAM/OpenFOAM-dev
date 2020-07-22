/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::AveragingMethods::Dual

Description
    Dual-mesh lagrangian averaging procedure.

    Point values are summed using the tetrahedral decomposition of the
    computational cells. Summation is done in the cells, and also in the
    tetrahedrons surrounding each point. The latter forms a type of dual mesh.
    The interpolation is weighted by proximity to the cell centre or point, as
    calculated by the barycentric coordinate within the tetrahedron.

    Values are interpolated linearly across the tetrahedron. Gradients are
    calculated directly from the point values using a first order finite
    element basis. The computed gradient is assumed constant over the
    tetrahedron.

SourceFiles
    Dual.C

\*---------------------------------------------------------------------------*/

#ifndef Dual_H
#define Dual_H

#include "AveragingMethod.H"
#include "pointMesh.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace AveragingMethods
{

/*---------------------------------------------------------------------------*\
                            Class Dual Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Dual
:
    public AveragingMethod<Type>
{
public:

    // Public Typedefs

        //- Gradient type
        typedef typename AveragingMethod<Type>::TypeGrad TypeGrad;


private:

    // Private data

        //- Volume of the cell-centered regions
        const Field<scalar>& volumeCell_;

        //- Volume of the point-centered regions
        Field<scalar> volumeDual_;

        //- Data on the cells
        Field<Type>& dataCell_;

        //- Data on the points
        Field<Type>& dataDual_;


    // Private static member functions

        //- Return the size of the FieldField parts
        static autoPtr<labelList> size(const fvMesh &mesh);


    // Private Member Functions

        //- Sync point data over processor boundaries
        void syncDualData();


public:

    //- Runtime type information
    TypeName("dual");


    //- Constructors

        //- Construct from components
        Dual
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        Dual(const Dual<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<AveragingMethod<Type>> clone() const
        {
            return autoPtr<AveragingMethod<Type>>
            (
                new Dual<Type>(*this)
            );
        }


    //- Destructor
    virtual ~Dual();


    //- Member Functions

        //- Add point value to interpolation
        void add
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const Type& value
        );

        //- Interpolate
        Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- Interpolate gradient
        TypeGrad interpolateGrad
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- Calculate the average
        void average();
        void average(const AveragingMethod<scalar>& weight);

        //- Return an internal field of the average
        tmp<Field<Type>> primitiveField() const;

        //- Return an internal field of the gradient
        tmp<Field<TypeGrad>> internalFieldGrad() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AveragingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Dual.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
