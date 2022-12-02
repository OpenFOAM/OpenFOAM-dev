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

\*---------------------------------------------------------------------------*/

#include "gaussLaplacianScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvLaplacianScheme(gaussLaplacianScheme)

#define declareFvmLaplacianScalarGamma(Type)                                   \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::fvMatrix<Foam::Type>>                                          \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvmLaplacian         \
(                                                                              \
    const SurfaceField<scalar>& gamma,                                         \
    const VolField<Type>& vf                                                   \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    SurfaceField<scalar> gammaMagSf                                            \
    (                                                                          \
        gamma*mesh.magSf()                                                     \
    );                                                                         \
                                                                               \
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected                         \
    (                                                                          \
        gammaMagSf,                                                            \
        this->tsnGradScheme_().deltaCoeffs(vf),                                \
        vf                                                                     \
    );                                                                         \
    fvMatrix<Type>& fvm = tfvm.ref();                                          \
                                                                               \
    if (this->tsnGradScheme_().corrected())                                    \
    {                                                                          \
        if (mesh.schemes().fluxRequired(vf.name()))                            \
        {                                                                      \
            fvm.faceFluxCorrectionPtr() = new                                  \
            SurfaceField<Type>                   \
            (                                                                  \
                gammaMagSf*this->tsnGradScheme_().correction(vf)               \
            );                                                                 \
                                                                               \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    *fvm.faceFluxCorrectionPtr()                               \
                )().primitiveField();                                          \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            fvm.source() -=                                                    \
                mesh.V()*                                                      \
                fvc::div                                                       \
                (                                                              \
                    gammaMagSf*this->tsnGradScheme_().correction(vf)           \
                )().primitiveField();                                          \
        }                                                                      \
    }                                                                          \
                                                                               \
    return tfvm;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<>                                                                     \
Foam::tmp<Foam::VolField<Foam::Type>> \
Foam::fv::gaussLaplacianScheme<Foam::Type, Foam::scalar>::fvcLaplacian         \
(                                                                              \
    const SurfaceField<scalar>& gamma,                                         \
    const VolField<Type>& vf                                                   \
)                                                                              \
{                                                                              \
    const fvMesh& mesh = this->mesh();                                         \
                                                                               \
    tmp<VolField<Type>> tLaplacian                                             \
    (                                                                          \
        fvc::div(gamma*this->tsnGradScheme_().snGrad(vf)*mesh.magSf())         \
    );                                                                         \
                                                                               \
    tLaplacian.ref().rename                                                    \
    (                                                                          \
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'                    \
    );                                                                         \
                                                                               \
    return tLaplacian;                                                         \
}


declareFvmLaplacianScalarGamma(scalar);
declareFvmLaplacianScalarGamma(vector);
declareFvmLaplacianScalarGamma(sphericalTensor);
declareFvmLaplacianScalarGamma(symmTensor);
declareFvmLaplacianScalarGamma(tensor);


// ************************************************************************* //
