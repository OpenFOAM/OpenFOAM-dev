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
    Foam::blended

Description
    linear/upwind blended interpolation scheme.

SourceFiles
    blended.C

\*---------------------------------------------------------------------------*/

#ifndef blended_H
#define blended_H

#include "limitedSurfaceInterpolationScheme.H"
#include "blendedSchemeBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blended
:
    public limitedSurfaceInterpolationScheme<Type>,
    public blendedSchemeBase<Type>
{
    // Private Data

        const scalar blendingFactor_;


public:

    //- Runtime type information
    TypeName("blended");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        blended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const scalar blendingFactor
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            blendingFactor_(blendingFactor)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        blended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is),
            blendingFactor_(readScalar(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        blended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            blendingFactor_(readScalar(is))
        {}

        //- Disallow default bitwise copy construction
        blended(const blended&) = delete;


    //- Destructor
    virtual ~blended()
    {}


    // Member Functions

        //- Return the face-based blending factor
        virtual tmp<surfaceScalarField> blendingFactor
        (
             const VolField<Type>& vf
        ) const
        {
            return surfaceScalarField::New
            (
                vf.name() + "BlendingFactor",
                this->mesh(),
                dimensionedScalar(dimless, blendingFactor_)
            );
        }

        //- Return the interpolation limiter
        virtual tmp<surfaceScalarField> limiter
        (
            const VolField<Type>&
        ) const
        {
            return surfaceScalarField::New
            (
                "blendedLimiter",
                this->mesh(),
                dimensionedScalar(dimless, 1 - blendingFactor_)
            );
        }

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const
        {
            return
                blendingFactor_*this->mesh().surfaceInterpolation::weights()
              + (1 - blendingFactor_)*pos0(this->faceFlux_);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blended&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
