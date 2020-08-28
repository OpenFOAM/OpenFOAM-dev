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
    Foam::multivariateScheme

Description
    Generic multi-variate discretisation scheme class which may be instantiated
    for any of the NVD, CNVD or NVDV schemes.

SourceFiles
    multivariateScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateScheme_H
#define multivariateScheme_H

#include "multivariateSurfaceInterpolationScheme.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class multivariateScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Scheme>
class multivariateScheme
:
    public multivariateSurfaceInterpolationScheme<Type>,
    public Scheme::LimiterType
{
    // Private Data

        const surfaceScalarField& faceFlux_;
        surfaceScalarField weights_;


public:

    //- Runtime type information
    TypeName("multivariateScheme");


    // Constructors

        //- Construct for field, faceFlux and Istream
        multivariateScheme
        (
            const fvMesh& mesh,
            const typename multivariateSurfaceInterpolationScheme<Type>::
                fieldTable& fields,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );

        //- Disallow default bitwise copy construction
        multivariateScheme(const multivariateScheme&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multivariateScheme&) = delete;

        //- surfaceInterpolationScheme sub-class returned by operator(field)
        //  which is used as the interpolation scheme for the field
        class fieldScheme
        :
            public multivariateSurfaceInterpolationScheme<Type>::
                fieldScheme
        {
            // Private Data

                const surfaceScalarField& weights_;

        public:

            // Constructors

                fieldScheme
                (
                    const GeometricField<Type, fvPatchField, volMesh>& field,
                    const surfaceScalarField& weights
                )
                :
                    multivariateSurfaceInterpolationScheme<Type>::
                        fieldScheme(field),
                    weights_(weights)
                {}


            // Member Functions

                //- Return the interpolation weighting factors
                tmp<surfaceScalarField> weights
                (
                    const GeometricField<Type, fvPatchField, volMesh>&
                ) const
                {
                    return weights_;
                }
        };

        tmp<surfaceInterpolationScheme<Type>> operator()
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const
        {
            return tmp<surfaceInterpolationScheme<Type>>
            (
                new fieldScheme(field, weights_)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLimitedMultivariateSurfaceInterpolationScheme(SS, LIMITER)         \
typedef multivariateScheme                                                     \
<                                                                              \
    scalar,                                                                    \
    LimitedScheme<scalar, LIMITER<NVDTVD>, limitFuncs::magSqr>                 \
>                                                                              \
    multivariateScheme##LIMITER##_;                                            \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        multivariateScheme##LIMITER##_,                                        \
        #SS,                                                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
multivariateSurfaceInterpolationScheme<scalar>::addIstreamConstructorToTable   \
<                                                                              \
    multivariateScheme                                                         \
    <                                                                          \
        scalar,                                                                \
        LimitedScheme<scalar, LIMITER<NVDTVD>, limitFuncs::magSqr>             \
    >                                                                          \
>                                                                              \
    addMultivariate##SS##ConstructorToTable_;


#define makeLLimitedMultivariateSurfaceInterpolationScheme\
(                                                                              \
    SS,                                                                        \
    LLIMITER,                                                                  \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC                                                                    \
)                                                                              \
typedef multivariateScheme                                                     \
<                                                                              \
    scalar,                                                                    \
    LimitedScheme<scalar, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>      \
>                                                                              \
    multivariateScheme##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;                 \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        multivariateScheme##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_,             \
        #SS,                                                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
multivariateSurfaceInterpolationScheme<scalar>::addIstreamConstructorToTable   \
<                                                                              \
    multivariateScheme                                                         \
    <                                                                          \
        scalar,                                                                \
        LimitedScheme<scalar, LLIMITER<LIMITER<NVDTVD>>, limitFuncs::LIMFUNC>  \
    >                                                                          \
>                                                                              \
    addMultivariate##SS##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multivariateScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
