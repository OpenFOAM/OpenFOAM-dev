/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::fv::EulerDdtScheme

Description
    Basic first-order Euler implicit/explicit ddt using only the current and
    previous time-step values.

SourceFiles
    EulerDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerDdtScheme_H
#define EulerDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class EulerDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerDdtScheme
:
    public ddtScheme<Type>
{
public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from mesh
        EulerDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        EulerDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        EulerDdtScheme(const EulerDdtScheme&) = delete;


    // Member Functions

        using ddtScheme<Type>::mesh;

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fvcDdt
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        virtual tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        virtual tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const EulerDdtScheme&) = delete;
};


template<>
tmp<surfaceScalarField> EulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> EulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> EulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> EulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EulerDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
