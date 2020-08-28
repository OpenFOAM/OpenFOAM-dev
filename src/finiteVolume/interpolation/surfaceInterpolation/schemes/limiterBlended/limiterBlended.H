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
    Foam::limiterBlended

Description
    Blends two specified schemes using the limiter function provided by a
    limitedSurfaceInterpolationScheme.

    The limited scheme is specified first followed by the scheme to be scaled
    by the limiter and then the scheme scaled by 1 - limiter e.g.

    div(phi,U)      Gauss limiterBlended vanLeer linear linearUpwind grad(U);

SourceFiles
    limiterBlended.C

\*---------------------------------------------------------------------------*/

#ifndef limiterBlended_H
#define limiterBlended_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class limiterBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limiterBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Limited scheme providing the limiter
        tmp<limitedSurfaceInterpolationScheme<Type>> tLimitedScheme_;

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;


public:

    //- Runtime type information
    TypeName("limiterBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        limiterBlended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tLimitedScheme_
            (
                limitedSurfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        limiterBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tLimitedScheme_
            (
                limitedSurfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        limiterBlended(const limiterBlended&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField blendingFactor
            (
                tLimitedScheme_().limiter(vf)
            );

            return
                blendingFactor*tScheme1_().weights(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().weights(vf);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate(const GeometricField<Type, fvPatchField, volMesh>& vf) const
        {
            surfaceScalarField blendingFactor
            (
                tLimitedScheme_().limiter(vf)
            );

            return
                blendingFactor*tScheme1_().interpolate(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().interpolate(vf);
        }


        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme1_().corrected() || tScheme2_().corrected();
        }


        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            surfaceScalarField blendingFactor
            (
                tLimitedScheme_().limiter(vf)
            );

            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        blendingFactor
                      * tScheme1_().correction(vf)
                      + (scalar(1) - blendingFactor)
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        blendingFactor
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return
                (
                    (scalar(1) - blendingFactor)
                  * tScheme2_().correction(vf)
                );
            }
            else
            {
                return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
                (
                    nullptr
                );
            }
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limiterBlended&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
