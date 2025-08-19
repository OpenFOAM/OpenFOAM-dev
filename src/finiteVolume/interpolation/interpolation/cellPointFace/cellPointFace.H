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
    Foam::interpolations::cellPointFace

Description
    Piecewise-linear interpolation method. Uses volPointInterpolation to create
    values on the points and linear surfaceInterpolation to create values at
    the face-centres. Then uses a face-centre tetrahedral decomposition to
    linearly interpolate between the values at the points, at the face-centres
    and at the cell-centres.

SourceFiles
    cellPointFace.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointFace_H
#define cellPointFace_H

#include "volPointInterpolation_interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolations
{

/*---------------------------------------------------------------------------*\
                 Class cellPointFace Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cellPointFace
:
    public fieldInterpolation<Type, cellPointFace<Type>>,
    public volPointInterpolation<Type>
{
    // Private Data

        //- Linearly interpolated volfield
        const SurfaceField<Type> psis_;


    // Private Member Functions

        //- Find the tetrahedron in which the position is. While searching for
        //  the tet, store the tet closest to the position. This way, if
        //  position is not inside any tet, choose the closest one.
        bool findTet
        (
            const vector& position,
            const label nFace,
            vector tetPoints[],
            label tetLabelCandidate[],
            label tetPointLabels[],
            scalar phi[],
            scalar phiCandidate[],
            label& closestFace,
            scalar& minDistance
        ) const;

        //- Find the triangle in which the position is, when the position is
        //  known to be on the face
        bool findTriangle
        (
            const vector& position,
            const label nFace,
            label tetPointLabels[],
            scalar phi[]
        ) const;


public:

    //- Runtime type information
    TypeName("cellPointFace");


    // Constructors

        //- Construct from components
        cellPointFace(const VolField<Type>& psi);

        //- Copy construct
        cellPointFace(const cellPointFace<Type>& i);

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new cellPointFace<Type>(*this)
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellPointFace.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
