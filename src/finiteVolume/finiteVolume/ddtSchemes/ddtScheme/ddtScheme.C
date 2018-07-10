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

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "surfaceInterpolate.H"
#include "fvMatrix.H"
#include "cyclicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<ddtScheme<Type>> ddtScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        InfoInFunction << "Constructing ddtScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Ddt scheme not specified" << endl << endl
            << "Valid ddt schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown ddt scheme " << schemeName << nl << nl
            << "Valid ddt schemes are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
ddtScheme<Type>::~ddtScheme()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> ddtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    NotImplemented;

    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        GeometricField<Type, fvPatchField, volMesh>::null()
    );
}


template<class Type>
tmp<fvMatrix<Type>> ddtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    NotImplemented;

    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()
            *vf.dimensions()*dimVol/dimTime
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> ddtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    NotImplemented;

    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::null()
    );
}



template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr
)
{
    tmp<surfaceScalarField> tddtCouplingCoeff = scalar(1)
      - min
        (
            mag(phiCorr)
           *mesh().time().deltaT()*mesh().deltaCoeffs()/mesh().magSf(),
            scalar(1)
        );

    surfaceScalarField& ddtCouplingCoeff = tddtCouplingCoeff.ref();

    surfaceScalarField::Boundary& ccbf =
        ddtCouplingCoeff.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
            U.boundaryField()[patchi].fixesValue()
         || isA<cyclicAMIFvPatch>(mesh().boundary()[patchi])
        )
        {
            ccbf[patchi] = 0.0;
        }
    }

    if (debug > 1)
    {
        InfoInFunction
            << "ddtCouplingCoeff mean max min = "
            << gAverage(ddtCouplingCoeff.primitiveField())
            << " " << gMax(ddtCouplingCoeff.primitiveField())
            << " " << gMin(ddtCouplingCoeff.primitiveField())
            << endl;
    }

    return tddtCouplingCoeff;
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const fluxFieldType& phiCorr,
    const volScalarField& rho
)
{
    return fvcDdtPhiCoeff(U, phi, phiCorr/fvc::interpolate(rho));
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return fvcDdtPhiCoeff(U, phi, phi - fvc::dotInterpolate(mesh().Sf(), U));
}


template<class Type>
tmp<surfaceScalarField> ddtScheme<Type>::fvcDdtPhiCoeff
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi,
    const volScalarField& rho
)
{
    return fvcDdtPhiCoeff
    (
        U,
        phi,
        (phi - fvc::dotInterpolate(mesh().Sf(), U))/fvc::interpolate(rho)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
