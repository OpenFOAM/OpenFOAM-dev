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

#include "slicedFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& completeField
)
:
    fvPatchField<Type>(p, iF, Field<Type>())
{
    // Set the fvPatchField to a slice of the given complete field
    UList<Type>::shallowCopy(p.patchSlice(completeField));
}


template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF, Field<Type>())
{}


template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false)
{
    NotImplemented;
}


template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{
    NotImplemented;
}


template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf.patch(), iF, Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::shallowCopy(ptf);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::slicedFvPatchField<Type>::clone() const
{
    return tmp<fvPatchField<Type>>
    (
        new slicedFvPatchField<Type>(*this)
    );
}


template<class Type>
Foam::slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>
    (
        ptf.patch(),
        ptf.internalField(),
        Field<Type>()
    )
{
    // Transfer the slice from the argument
    UList<Type>::shallowCopy(ptf);
}


template<class Type>
Foam::tmp<Foam::fvPatchField<Type>>
Foam::slicedFvPatchField<Type>::clone
(
    const DimensionedField<Type, volMesh>& iF
) const
{
    return tmp<fvPatchField<Type>>
    (
        new slicedFvPatchField<Type>(*this, iF)
    );
}


template<class Type>
Foam::slicedFvPatchField<Type>::~slicedFvPatchField()
{
    // Set the fvPatchField storage pointer to nullptr before its destruction
    // to protect the field it a slice of.
    UList<Type>::shallowCopy(UList<Type>(nullptr, 0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::slicedFvPatchField<Type>::snGrad() const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
void Foam::slicedFvPatchField<Type>::updateCoeffs()
{
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::patchInternalField() const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
void Foam::slicedFvPatchField<Type>::patchInternalField(Field<Type>&) const
{
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::patchNeighbourField() const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::gradientInternalCoeffs() const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    NotImplemented;

    return Field<Type>::null();
}


template<class Type>
void Foam::slicedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
