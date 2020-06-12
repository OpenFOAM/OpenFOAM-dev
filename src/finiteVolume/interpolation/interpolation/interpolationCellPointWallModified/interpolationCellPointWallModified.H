/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::interpolationCellPointWallModified

Description
    As interpolationCellPoint, but with the point field modified on wall faces.

    This method is defined only for vectors. The point field is extrapolated
    from the cells to the wall faces, and then rotated towards the (reverse)
    point normal so that the vectors do not point out of the domain. The result
    is also scaled so so if the necessary rotation exceeds 90 degrees, it is
    clamped to zero.

    This prevents unresolvable drag-rebound couplings when applied to the
    velocity interpolation in a Lagrangian simulation.

\*---------------------------------------------------------------------------*/

#ifndef interpolationCellPointWallModified_H
#define interpolationCellPointWallModified_H

#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class interpolationCellPointWallModified Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationCellPointWallModified
:
    public interpolationCellPoint<Type>
{
    // Private Member Functions

        //- Compute the point field for a type. This just throws an error. This
        //  interpolation method is only defined for vector fields.
        template <class TYPE>
        tmp<PointField<TYPE>> calcPointField
        (
            const VolField<TYPE>& psi
        ) const;

        //- Compute the point field for a vector type
        tmp<pointVectorField> calcPointField(const volVectorField& psi) const;


public:

    //- Runtime type information
    TypeName("cellPointWallModified");


    // Constructors

        //- Construct from components
        interpolationCellPointWallModified
        (
            const VolField<Type>& psi
        );

        //- Copy construct
        interpolationCellPointWallModified
        (
            const interpolationCellPointWallModified<Type>& i
        );

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new interpolationCellPointWallModified<Type>(*this)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationCellPointWallModified.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
