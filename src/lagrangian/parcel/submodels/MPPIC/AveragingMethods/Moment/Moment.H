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
    Foam::AveragingMethods::Moment

Description
    Moment lagrangian averaging procedure.

    Point values and moments from the cell centroid are summed over
    computational cells. A linear function is generated which has the same
    integrated moment as that of the point data.

    The computed linear function is used to interpolate values within a cell.
    The gradient is calculated from the coefficients of the function, and is
    assumed constant over the cell.

SourceFiles
    Moment.C

\*---------------------------------------------------------------------------*/

#ifndef Moment_H
#define Moment_H

#include "AveragingMethod.H"
#include "pointMesh.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace AveragingMethods
{

/*---------------------------------------------------------------------------*\
                           Class Moment Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Moment
:
    public AveragingMethod<Type>
{
public:

    // Public Typedefs

        //- Gradient type
        typedef typename AveragingMethod<Type>::TypeGrad TypeGrad;


private:

    // Private data

        //- Data mean
        Field<Type>& data_;

        //- X-data moment
        Field<Type>& dataX_;

        //- Y-data moment
        Field<Type>& dataY_;

        //- Z-data moment
        Field<Type>& dataZ_;

        //- Transform tensor from moment to gradient
        Field<symmTensor> transform_;

        //- Length scale for moment values
        Field<scalar> scale_;


    // Private Member Functions

        //- Re-calculate gradient
        virtual void updateGrad();


public:

    //- Runtime type information
    TypeName("moment");


    //- Constructors

        //- Construct from components
        Moment
        (
            const IOobject& io,
            const dictionary& dict,
            const fvMesh &mesh
        );

        //- Construct a copy
        Moment(const Moment<Type>& am);

        //- Construct and return a clone
        virtual autoPtr<AveragingMethod<Type>> clone() const
        {
            return autoPtr<AveragingMethod<Type>>
            (
                new Moment<Type>(*this)
            );
        }


    //- Destructor
    virtual ~Moment();


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
    #include "Moment.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
