/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::interpolationCellPoint

Description
    Given cell centre values and point (vertex) values decompose into
    tetrahedra and linear interpolate within them.

\*---------------------------------------------------------------------------*/

#ifndef interpolationCellPoint_H
#define interpolationCellPoint_H

#include "interpolation.H"
#include "cellPointWeight.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class interpolationCellPoint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationCellPoint
:
    public fieldInterpolation<Type, interpolationCellPoint<Type>>
{
protected:

    // Protected data

        //- Interpolated volfield
        const GeometricField<Type, pointPatchField, pointMesh> psip_;


public:

    //- Runtime type information
    TypeName("cellPoint");


    // Constructors

        //- Construct from components
        interpolationCellPoint
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        //- Construct from components
        interpolationCellPoint
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            tmp<GeometricField<Type, pointPatchField, pointMesh>> psip
        );


    // Member Functions

        //- Interpolate field for the given cellPointWeight
        inline Type interpolate(const cellPointWeight& cpw) const;

        //- Interpolate field to the given point in the given cell
        inline Type interpolate
        (
            const vector& position,
            const label celli,
            const label facei = -1
        ) const;

        //- Interpolate field to the given coordinates in the tetrahedron
        //  defined by the given indices.
        inline Type interpolate
        (
            const barycentric& coordinates,
            const tetIndices& tetIs,
            const label facei = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationCellPointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationCellPoint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
