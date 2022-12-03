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
    Foam::fv::gaussLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussLaplacianScheme_H
#define gaussLaplacianScheme_H

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GType>
class gaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{
    // Private Member Functions

        tmp<SurfaceField<Type>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const VolField<Type>&
        );

        //- Disallow default bitwise copy construction
        gaussLaplacianScheme(const gaussLaplacianScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const gaussLaplacianScheme&) = delete;


public:

    //- Runtime type information
    TypeName("Gauss");


    // Constructors

        //- Construct null
        gaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is)
        {}

        //- Construct from mesh, interpolation and snGradScheme schemes
        gaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}


    //- Destructor
    virtual ~gaussLaplacianScheme()
    {}


    // Member Functions

        static tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& deltaCoeffs,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcLaplacian
        (
            const VolField<Type>&
        );

        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );

        tmp<VolField<Type>> fvcLaplacian
        (
            const SurfaceField<GType>&,
            const VolField<Type>&
        );
};


// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma

#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>> gaussLaplacianScheme<Type, scalar>::fvmLaplacian           \
(                                                                              \
    const SurfaceField<scalar>&,                 \
    const VolField<Type>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<VolField<Type>>                               \
gaussLaplacianScheme<Type, scalar>::fvcLaplacian                               \
(                                                                              \
    const SurfaceField<scalar>&,                 \
    const VolField<Type>&                         \
);


defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
