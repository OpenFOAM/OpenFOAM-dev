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

#include "steadyStateDdtScheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        "ddt("+dt.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            dt.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        "ddt("+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            vf.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            rho.dimensions()*vf.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            rho.dimensions()*vf.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            rho.dimensions()*vf.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        mesh(),
        dimensioned<typename flux<Type>::type>
        (
            "0",
            Uf.dimensions()*dimArea/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        mesh(),
        dimensioned<typename flux<Type>::type>
        (
            "0",
            phi.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return fluxFieldType::New
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
        mesh(),
        dimensioned<typename flux<Type>::type>
        (
            "0",
            Uf.dimensions()*dimArea/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return fluxFieldType::New
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh(),
        dimensioned<typename flux<Type>::type>
        (
            "0",
            phi.dimensions()/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<surfaceScalarField> steadyStateDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return surfaceScalarField::New
    (
        "meshPhi",
        mesh(),
        dimensionedScalar(dimVolume/dimTime, 0)
    );
}


template<class Type>
tmp<scalarField> steadyStateDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const label patchi
)
{
    return tmp<scalarField>
    (
        new scalarField(mesh().boundary()[patchi].size(), 0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
