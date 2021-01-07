/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type, class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
Foam::dimensionSet Foam::fv::option::sourceDims
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
Foam::dimensionSet Foam::fv::option::sourceDims
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const dimensionSet& ds
)
{
    return field.dimensions()*ds;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, class ... AlphaRhoFieldTypes>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::source
(
    GeometricField<Type, fvPatchField, volMesh>& field,
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

    const label fieldi = applyToField(fieldName);

    if (fieldi != -1)
    {
        addSup(alphaRhos ..., mtx, fieldi);
    }

    return tmtx;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::option::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
) const
{
    return source(field, fieldName, dimVolume/dimTime);
}


// ************************************************************************* //
