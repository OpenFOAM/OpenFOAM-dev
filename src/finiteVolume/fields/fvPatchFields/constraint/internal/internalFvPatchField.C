/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "internalFvPatchField.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::internalFvPatchField<Type>::internalFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::internalFvPatchField<Type>::internalFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false)
{
    if (!isType<internalFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    fvPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
Foam::internalFvPatchField<Type>::internalFvPatchField
(
    const internalFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<internalFvPatch>(p))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}


template<class Type>
Foam::internalFvPatchField<Type>::internalFvPatchField
(
    const internalFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::internalFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::operator==(this->patchInternalField());
    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    if (this->patch().size())
    {
        FatalErrorInFunction
            << "attempt to create matrix coefficients for field "
            << this->internalField().name()
            << " on non-empty '" << typeName << "' patch "
            << this->patch().name()
            << exit(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    if (this->patch().size())
    {
        FatalErrorInFunction
            << "attempt to create matrix coefficients for field "
            << this->internalField().name()
            << " on non-empty '" << typeName << "' patch "
            << this->patch().name()
            << exit(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalFvPatchField<Type>::gradientInternalCoeffs() const
{
    if (this->patch().size())
    {
        FatalErrorInFunction
            << "attempt to create matrix coefficients for field "
            << this->internalField().name()
            << " on non-empty '" << typeName << "' patch "
            << this->patch().name()
            << exit(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::internalFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    if (this->patch().size())
    {
        FatalErrorInFunction
            << "attempt to create matrix coefficients for field "
            << this->internalField().name()
            << " on non-empty '" << typeName << "' patch "
            << this->patch().name()
            << exit(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


// ************************************************************************* //
