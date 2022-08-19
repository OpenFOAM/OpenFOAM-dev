/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type, class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
Foam::dimensionSet Foam::fvModel::sourceDims
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const dimensionSet& ds,
    const AlphaRhoFieldType& alphaRho,
    const AlphaRhoFieldTypes& ... alphaRhos
)
{
    return alphaRho.dimensions()*sourceDims(field, ds, alphaRhos ...);
}


template<class Type>
Foam::dimensionSet Foam::fvModel::sourceDims
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const dimensionSet& ds
)
{
    return field.dimensions()*ds;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class ... AlphaRhoFieldTypes>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName,
    const dimensionSet& ds,
    const AlphaRhoFieldTypes& ... alphaRhos
) const
{
    tmp<fvMatrix<Type>> tmtx
    (
        new fvMatrix<Type>
        (
            field,
            sourceDims(field, ds, alphaRhos ...)
        )
    );
    fvMatrix<Type>& mtx = tmtx.ref();

    if (addsSupToField(fieldName))
    {
        addSup(alphaRhos ..., mtx, fieldName);
    }

    return tmtx;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->source(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->source(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime, rho);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->source(alpha, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime, alpha, rho);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->source(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    volScalarField one
    (
        IOobject
        (
            "one",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    );

    return this->source(alpha, one, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::source
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->source(rho, field, field.name());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::d2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->d2dt2(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fvModel::d2dt2
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/sqr(dimTime));
}

// ************************************************************************* //
