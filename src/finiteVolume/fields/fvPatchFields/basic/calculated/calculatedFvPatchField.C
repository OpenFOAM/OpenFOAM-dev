/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "calculatedFvPatchField.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    Type calculatedFvPatchFieldNaN()
    {
        return pTraits<Type>::nan;
    }

    template<>
    inline label calculatedFvPatchFieldNaN<label>()
    {
        return -labelMax;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::fvPatchField<Type>::calculatedType()
{
    return calculatedFvPatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<Type>(p, iF, dict, valueRequired)
{}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, false)
{
    if (mappingRequired)
    {
        mapper(*this, ptf, calculatedFvPatchFieldNaN<Type>());
    }
}


template<class Type>
Foam::calculatedFvPatchField<Type>::calculatedFvPatchField
(
    const calculatedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::fvPatchField<Type>::NewCalculatedType
(
    const fvPatch& p
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()
        (
            p,
            DimensionedField<Type, volMesh>::null()
        );
    }
    else
    {
        return tmp<fvPatchField<Type>>
        (
            new calculatedFvPatchField<Type>
            (
                p,
                DimensionedField<Type, volMesh>::null()
            )
        );
    }
}


template<class Type>
template<class Type2>
Foam::tmp<Foam::fvPatchField<Type>> Foam::fvPatchField<Type>::NewCalculatedType
(
    const fvPatchField<Type2>& pf
)
{
    return NewCalculatedType(pf.patch());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::calculatedFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    mapper(*this, ptf, calculatedFvPatchFieldNaN<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorInFunction
        << "cannot be called for a calculatedFvPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->internalField().name()
        << " in file " << this->internalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << abort(FatalError);

    return *this;
}


template<class Type>
void Foam::calculatedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
