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
    Foam::interpolationCellPointFace

Description
    Foam::interpolationCellPointFace

\*---------------------------------------------------------------------------*/

#ifndef interpolationCellPointFace_H
#define interpolationCellPointFace_H

#include "interpolationVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class interpolationCellPointFace Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interpolationCellPointFace
:
    public fieldInterpolation<Type, interpolationCellPointFace<Type>>,
    public interpolationVolPointInterpolation<Type>
{
    // Private Data

        //- Linearly interpolated volfield
        const SurfaceField<Type> psis_;

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
        interpolationCellPointFace
        (
            const VolField<Type>& psi
        );

        //- Copy construct
        interpolationCellPointFace
        (
            const interpolationCellPointFace<Type>& i
        );

        //- Clone
        virtual autoPtr<interpolation<Type>> clone() const
        {
            return autoPtr<interpolation<Type>>
            (
                new interpolationCellPointFace<Type>(*this)
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "interpolationCellPointFace.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
