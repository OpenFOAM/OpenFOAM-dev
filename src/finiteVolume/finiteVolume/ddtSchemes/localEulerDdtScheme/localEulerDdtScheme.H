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
    Foam::fv::localEulerDdtScheme

Description
    Local time-step first-order Euler implicit/explicit ddt.

    The reciprocal of the local time-step field is looked-up from the database.

    This scheme should only be used for steady-state computations using
    transient codes where local time-stepping is preferably to under-relaxation
    for transport consistency reasons.

See also
    Foam::fv::CoEulerDdtScheme

SourceFiles
    localEulerDdt.C
    localEulerDdtScheme.C
    localEulerDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef localEulerDdtScheme_H
#define localEulerDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class localEulerDdt Declaration
\*---------------------------------------------------------------------------*/

class localEulerDdt
{
public:

    //- Name of the reciprocal local time-step field
    static word rDeltaTName;

    //- Name of the reciprocal local face time-step field
    static word rDeltaTfName;

    //- Name of the reciprocal local sub-cycling time-step field
    static word rSubDeltaTName;


    // Constructors

        localEulerDdt()
        {}


    // Member Functions

        //- Return true if LTS is enabled
        static bool enabled(const fvMesh& mesh);

        //- Return the reciprocal of the local time-step
        //  looked-up from the objectRegistry
        static const volScalarField& localRDeltaT(const fvMesh& mesh);

        //- Return the reciprocal of the local face time-step
        //  looked-up from the objectRegistry
        static const surfaceScalarField& localRDeltaTf(const fvMesh& mesh);

        //- Calculate and return the reciprocal of the local sub-cycling
        //  time-step
        static tmp<volScalarField> localRSubDeltaT
        (
            const fvMesh& mesh,
            const label nAlphaSubCycles
        );
};


/*---------------------------------------------------------------------------*\
                     Class localEulerDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localEulerDdtScheme
:
    public localEulerDdt,
    public fv::ddtScheme<Type>
{
    // Private Member Functions

        //- Return the reciprocal of the local time-step
        const volScalarField& localRDeltaT() const;

        //- Return the reciprocal of the local face time-step
        const surfaceScalarField& localRDeltaTf() const;


public:

    //- Runtime type information
    TypeName("localEuler");


    // Constructors

        //- Construct from mesh
        localEulerDdtScheme(const fvMesh& mesh)
        :
            ddtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        localEulerDdtScheme(const fvMesh& mesh, Istream& is)
        :
            ddtScheme<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        localEulerDdtScheme(const localEulerDdtScheme&) = delete;


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

        // using ddtScheme<Type>::fvcDdtPhiCoeff;

        // virtual tmp<surfaceScalarField> fvcDdtPhiCoeff
        // (
        //     const GeometricField<Type, fvPatchField, volMesh>& U,
        //     const fluxFieldType& phi,
        //     const fluxFieldType& phiCorr
        // );

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
        void operator=(const localEulerDdtScheme&) = delete;
};


template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> localEulerDdtScheme<scalar>::fvcDdtPhiCorr
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
    #include "localEulerDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
