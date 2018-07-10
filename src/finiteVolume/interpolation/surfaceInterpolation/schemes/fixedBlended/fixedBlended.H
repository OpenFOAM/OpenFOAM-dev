/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fixedBlended

Description
    Two-scheme fixed-blending differencing scheme.

    Similar to localBlended but uses a single (global) constant blending
    factor. The factor applies to the first scheme and 1-factor to the
    second scheme.

Note
    Although a blending factor of 0 and 1 is permitted, it is more efficient
    just to use the underlying scheme directly.

SourceFiles
    fixedBlended.C

\*---------------------------------------------------------------------------*/

#ifndef fixedBlended_H
#define fixedBlended_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fixedBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fixedBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private data

        const scalar blendingFactor_;

    // Private Member Functions

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type>> tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type>> tScheme2_;


        //- Disallow default bitwise copy construct
        fixedBlended(const fixedBlended&);

        //- Disallow default bitwise assignment
        void operator=(const fixedBlended&);


public:

    //- Runtime type information
    TypeName("fixedBlended");


    // Constructors

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        fixedBlended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            blendingFactor_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {
            if (blendingFactor_ < 0 || blendingFactor_ > 1)
            {
                FatalIOErrorInFunction(is)
                    << "coefficient = " << blendingFactor_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
            if (surfaceInterpolationScheme<Type>::debug)
            {
                Info<<"fixedBlended: " << blendingFactor_
                    << "*" << tScheme1_().type()
                    << " + (1-" << blendingFactor_ << ")*"
                    << tScheme2_().type()
                    <<endl;
            }
        }


        //- Construct from mesh, faceFlux and Istream
        fixedBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            blendingFactor_(readScalar(is)),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {
            if (blendingFactor_ < 0 || blendingFactor_ > 1)
            {
                FatalIOErrorInFunction(is)
                    << "coefficient = " << blendingFactor_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
            if (surfaceInterpolationScheme<Type>::debug)
            {
                Info<<"fixedBlended: " << blendingFactor_
                    << "*" << tScheme1_().type()
                    << " + (1-" << blendingFactor_ << ")*"
                    << tScheme2_().type()
                    <<endl;
            }
        }


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField>
        weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return
                blendingFactor_*tScheme1_().weights(vf)
              + (scalar(1) - blendingFactor_)*tScheme2_().weights(vf);
        }


        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return
                blendingFactor_*tScheme1_().interpolate(vf)
              + (scalar(1) - blendingFactor_)*tScheme2_().interpolate(vf);
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
            if (tScheme1_().corrected())
            {
                if (tScheme2_().corrected())
                {
                    return
                    (
                        blendingFactor_
                      * tScheme1_().correction(vf)
                      + (scalar(1) - blendingFactor_)
                      * tScheme2_().correction(vf)
                    );
                }
                else
                {
                    return
                    (
                        blendingFactor_
                      * tScheme1_().correction(vf)
                    );
                }
            }
            else if (tScheme2_().corrected())
            {
                return
                (
                    (scalar(1) - blendingFactor_)
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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
