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
    Foam::outletStabilised

Description
    Outlet-stabilised interpolation scheme which applies upwind interpolation
    to the faces of the cells adjacent to outlets.

    This is particularly useful to stabilise the velocity at entrainment
    boundaries for LES cases using linear or other centred interpolation
    schemes.

SourceFiles
    outletStabilised.C

\*---------------------------------------------------------------------------*/

#ifndef outletStabilised_H
#define outletStabilised_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class outletStabilised Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletStabilised
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const surfaceScalarField& faceFlux_;
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


public:

    //- Runtime type information
    TypeName("outletStabilised");


    // Constructors

        //- Construct from mesh and Istream
        outletStabilised
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            ),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux_, is)
            )
        {}


        //- Construct from mesh, faceFlux and Istream
        outletStabilised
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        outletStabilised(const outletStabilised&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        inline tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme_().corrected();
        }

        //- Return the explicit correction to the face-interpolate
        //  set to zero on the near-boundary faces where upwind is applied
        inline virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const outletStabilised&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
