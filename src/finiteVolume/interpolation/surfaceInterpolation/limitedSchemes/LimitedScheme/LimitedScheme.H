/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::LimitedScheme

Description
    Class to create NVD/TVD limited weighting-factors.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied the central-differencing weighting factor, the face-flux, the
    cell and face gradients (from which the normalised variable
    distribution may be created) and the cell centre distance.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    LimitedScheme.C

\*---------------------------------------------------------------------------*/

#ifndef LimitedScheme_H
#define LimitedScheme_H

#include "limitedSurfaceInterpolationScheme.H"
#include "LimitFuncs.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class LimitedScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Limiter, template<class> class LimitFunc>
class LimitedScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public Limiter
{
    // Private Member Functions

        //- Calculate the limiter
        void calcLimiter
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            surfaceScalarField& limiterField
        ) const;

        //- Disallow default bitwise copy construct
        LimitedScheme(const LimitedScheme&);

        //- Disallow default bitwise assignment
        void operator=(const LimitedScheme&);


public:

    //- Runtime type information
    TypeName("LimitedScheme");

    typedef Limiter LimiterType;

    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        LimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const Limiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(weight)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        LimitedScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is),
            Limiter(is)
        {}

        //- Construct from mesh, faceFlux and Istream
        LimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(is)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeLimitedSurfaceInterpolationTypeScheme\
(                                                                             \
    SS,                                                                       \
    LIMITER,                                                                  \
    NVDTVD,                                                                   \
    LIMFUNC,                                                                  \
    TYPE                                                                      \
)                                                                             \
                                                                              \
typedef LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC>             \
    LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_;                         \
defineTemplateTypeNameAndDebugWithName                                        \
    (LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);               \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >                  \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >                  \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                      \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable            \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >                  \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                   \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable        \
<LimitedScheme<TYPE, LIMITER<NVDTVD>, limitFuncs::LIMFUNC> >                  \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


#define makeLimitedSurfaceInterpolationScheme(SS, LIMITER)                    \
                                                                              \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,scalar)    \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,vector)    \
makeLimitedSurfaceInterpolationTypeScheme                                     \
(                                                                             \
    SS,                                                                       \
    LIMITER,                                                                  \
    NVDTVD,                                                                   \
    magSqr,                                                                   \
    sphericalTensor                                                           \
)                                                                             \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,symmTensor)\
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDTVD,magSqr,tensor)


#define makeLimitedVSurfaceInterpolationScheme(SS, LIMITER)                   \
makeLimitedSurfaceInterpolationTypeScheme(SS,LIMITER,NVDVTVDV,null,vector)


#define makeLLimitedSurfaceInterpolationTypeScheme\
(                                                                             \
    SS,                                                                       \
    LLIMITER,                                                                 \
    LIMITER,                                                                  \
    NVDTVD,                                                                   \
    LIMFUNC,                                                                  \
    TYPE                                                                      \
)                                                                             \
                                                                              \
typedef LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC>  \
    LimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;               \
defineTemplateTypeNameAndDebugWithName                                        \
    (LimitedScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);     \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                   \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> >       \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                          \
                                                                              \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable               \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> >       \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                      \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable            \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> >       \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                   \
                                                                              \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable        \
<LimitedScheme<TYPE, LLIMITER<LIMITER<NVDTVD> >, limitFuncs::LIMFUNC> >       \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LimitedScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
