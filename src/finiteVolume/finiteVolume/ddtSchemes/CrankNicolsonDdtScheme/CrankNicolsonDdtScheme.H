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
    Foam::fv::CrankNicolsonDdtScheme

Description
    Second-oder Crank-Nicolson implicit ddt using the current and
    previous time-step fields as well as the previous time-step ddt.

    The Crank-Nicolson scheme is often unstable for complex flows in complex
    geometries and it is necessary to "off-centre" the scheme to stabilize it
    while retaining greater temporal accuracy than the first-order
    Euler-implicit scheme.  Off-centering is specified via the mandatory
    coefficient \c ocCoeff in the range [0,1] following the scheme name e.g.
    \verbatim
    ddtSchemes
    {
        default         CrankNicolson 0.9;
    }
    \endverbatim
    or with an optional "ramp" function to transition from the Euler scheme to
    Crank-Nicolson over a initial period to avoid start-up problems, e.g.
    \verbatim
    ddtSchemes
    {
        default         CrankNicolson
        ocCoeff
        {
            type scale;
            scale linearRamp;
            duration 0.01;
            value 0.9;
        };
    }
    \endverbatim
    With a coefficient of 1 the scheme is fully centred and second-order,
    with a coefficient of 0 the scheme is equivalent to Euler-implicit.
    A coefficient of 0.9 has been found to be suitable for a range of cases for
    which higher-order accuracy in time is needed and provides similar accuracy
    and stability to the backward scheme.  However, the backward scheme has
    been found to be more robust and provides formal second-order accuracy in
    time.

    The advantage of the Crank-Nicolson scheme over backward is that only the
    new and old-time values are needed, the additional terms relating to the
    fluxes and sources are evaluated at the mid-point of the time-step which
    provides the opportunity to limit the fluxes in such a way as to ensure
    boundedness while maintaining greater accuracy in time compared to the
    Euler-implicit scheme.  This approach is now used with MULES in the
    interFoam family of solvers.  Boundedness cannot be guaranteed with the
    backward scheme.

Note
    The Crank-Nicolson coefficient for the implicit part of the RHS is related
    to the off-centering coefficient by
    \verbatim
        cnCoeff = 1.0/(1.0 + ocCoeff);
    \endverbatim

See also
    Foam::Euler
    Foam::backward

SourceFiles
    CrankNicolsonDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolsonDdtScheme_H
#define CrankNicolsonDdtScheme_H

#include "ddtScheme.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class CrankNicolsonDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CrankNicolsonDdtScheme
:
    public fv::ddtScheme<Type>
{
    // Private Data

        //- Class to store the ddt0 fields on the objectRegistry for use in the
        //  next time-step.  The start-time index of the CN scheme is also
        //  stored to help handle the transition from Euler to CN
        template<class GeoField>
        class DDt0Field
        :
            public GeoField
        {
            label startTimeIndex_;

        public:

            //- Constructor from file for restart.
            DDt0Field
            (
                const IOobject& io,
                const fvMesh& mesh
            );

            //- Constructor from components, initisalised to zero with given
            //  dimensions.
            DDt0Field
            (
                const IOobject& io,
                const fvMesh& mesh,
                const dimensioned<typename GeoField::value_type>& dimType
            );

            //- Return the start-time index
            label startTimeIndex() const;

            //- Cast to the underlying GeoField
            GeoField& operator()();

            //- Assignment to a GeoField
            void operator=(const GeoField& gf);
        };


        //- Off-centering coefficient function
        //  1 -> CN, less than one blends with EI
        autoPtr<Function1<scalar>> ocCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CrankNicolsonDdtScheme(const CrankNicolsonDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const CrankNicolsonDdtScheme&);

        template<class GeoField>
        DDt0Field<GeoField>& ddt0_
        (
            const word& name,
            const dimensionSet& dims
        );

        //- Check if the ddt0 needs to be evaluated for this time-step
        template<class GeoField>
        bool evaluate(const DDt0Field<GeoField>& ddt0) const;

        //- Return the coefficient for Euler scheme for the first time-step
        //  for and CN thereafter
        template<class GeoField>
        scalar coef_(const DDt0Field<GeoField>&) const;

        //- Return the old time-step coefficient for Euler scheme for the
        //  second time-step and for CN thereafter
        template<class GeoField>
        scalar coef0_(const DDt0Field<GeoField>&) const;

        //- Return the reciprocal time-step coefficient for Euler for the
        //  first time-step and CN thereafter
        template<class GeoField>
        dimensionedScalar rDtCoef_(const DDt0Field<GeoField>&) const;

        //- Return the reciprocal old time-step coefficient for Euler for the
        //  second time-step and CN thereafter
        template<class GeoField>
        dimensionedScalar rDtCoef0_(const DDt0Field<GeoField>&) const;

        //- Return ddt0 multiplied by the off-centreing coefficient
        template<class GeoField>
        tmp<GeoField> offCentre_(const GeoField& ddt0) const;


public:

    //- Runtime type information
    TypeName("CrankNicolson");


    // Constructors

        //- Construct from mesh
        CrankNicolsonDdtScheme(const fvMesh& mesh);

        //- Construct from mesh and Istream
        CrankNicolsonDdtScheme(const fvMesh& mesh, Istream& is);


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return fv::ddtScheme<Type>::mesh();
        }

        //- Return the current off-centreing coefficient
        scalar ocCoeff() const
        {
            return ocCoeff_->value(mesh().time().value());
        }

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensioned<Type>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );

        tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
        );

        tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi
        );


        tmp<surfaceScalarField> meshPhi
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};


template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtPhiCorr
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
    #include "CrankNicolsonDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
