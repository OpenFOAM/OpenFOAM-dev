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
    Foam::interpolations::cell

Description
    Basic piecewise-constant interpolation method. Interpolates to a point by
    using the value in the cell containing the point.

SourceFiles
    cell_interpolation.C

\*---------------------------------------------------------------------------*/

#ifndef cell_interpolation_H
#define cell_interpolation_H

#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolations
{

/*---------------------------------------------------------------------------*\
                            Class cell Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cell
:
    public fieldInterpolation<Type, cell<Type>>
{
public:

    //- Runtime type information
    TypeName("cell");


    // Constructors

        //- Construct from components
        cell(const VolField<Type>& psi);

        //- Copy construct
        cell(const cell<Type>& i);

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new cell<Type>(*this)
            );
        }


    // Member Functions

        //- Inherit interpolate from interpolation
        using interpolation<Type>::interpolate;

        //- Interpolate field to the given point in the given cell
        Type interpolate
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;

        //- Interpolate field to the given coordinates in the tetrahedron
        //  defined by the given indices. This is an optimisation which skips
        //  calculating the position, as cell interpolation doesn't need it.
        inline Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const
        {
            return interpolate(vector::zero, tetIs.cell(), facei);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cell_interpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
