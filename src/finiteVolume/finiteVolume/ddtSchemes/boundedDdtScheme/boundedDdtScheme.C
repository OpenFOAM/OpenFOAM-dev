/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "boundedDdtScheme.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
boundedDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return scheme_().fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
boundedDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
boundedDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
boundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvcDdt(rho, vf) - fvc::ddt(rho)*vf;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
boundedDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvcDdt(alpha, rho, vf) - fvc::ddt(alpha, rho)*vf;
}


template<class Type>
tmp<fvMatrix<Type> >
boundedDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvmDdt(vf);
}


template<class Type>
tmp<fvMatrix<Type> >
boundedDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvmDdt(rho, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
boundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvmDdt(rho, vf) - fvm::Sp(fvc::ddt(rho), vf);
}


template<class Type>
tmp<fvMatrix<Type> >
boundedDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().fvmDdt(alpha, rho, vf) - fvm::Sp(fvc::ddt(alpha, rho), vf);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return scheme_().fvcDdtUfCorr(U, Uf);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return scheme_().fvcDdtPhiCorr(U, phi);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return scheme_().fvcDdtUfCorr(rho, U, Uf);
}


template<class Type>
tmp<typename boundedDdtScheme<Type>::fluxFieldType>
boundedDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return scheme_().fvcDdtPhiCorr(rho, U, phi);
}


template<class Type>
tmp<surfaceScalarField> boundedDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_().meshPhi(vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
