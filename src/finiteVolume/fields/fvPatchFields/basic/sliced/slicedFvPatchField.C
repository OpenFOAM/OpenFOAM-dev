/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& completeField
)
:
    fvPatchField<Type>(p, iF, Field<Type>())
{
    // Set the fvPatchField to a slice of the given complete field
    UList<Type>::operator=(p.patchSlice(completeField));
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF, Field<Type>())
{}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "slicedFvPatchField(const slicedFvPatchField<Type>&, "
        "const fvPatch&, const Field<Type>&, const fvPatchFieldMapper&)"
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict)
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "slicedFvPatchField(const Field<Type>&, const dictionary&)"
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf.patch(), iF, Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::operator=(ptf);
}

template<class Type>
tmp<fvPatchField<Type> > slicedFvPatchField<Type>::clone() const
{
    return tmp<fvPatchField<Type> >
    (
        new slicedFvPatchField<Type>(*this)
    );
}


template<class Type>
slicedFvPatchField<Type>::slicedFvPatchField
(
    const slicedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>
    (
        ptf.patch(),
        ptf.dimensionedInternalField(),
        Field<Type>()
    )
{
    // Transfer the slice from the argument
    UList<Type>::operator=(ptf);
}


template<class Type>
tmp<fvPatchField<Type> > slicedFvPatchField<Type>::clone
(
    const DimensionedField<Type, volMesh>& iF
) const
{
    return tmp<fvPatchField<Type> >
    (
        new slicedFvPatchField<Type>(*this, iF)
    );
}


template<class Type>
slicedFvPatchField<Type>::~slicedFvPatchField<Type>()
{
    // Set the fvPatchField storage pointer to NULL before its destruction
    // to protect the field it a slice of.
    UList<Type>::operator=(UList<Type>(NULL, 0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::snGrad() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "snGrad()"
    );

    return Field<Type>::null();
}


template<class Type>
void slicedFvPatchField<Type>::updateCoeffs()
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "updateCoeffs()"
    );
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchInternalField() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchInternalField()"
    );

    return Field<Type>::null();
}


template<class Type>
void slicedFvPatchField<Type>::patchInternalField(Field<Type>&) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchInternalField(Field<Type>&)"
    );
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchNeighbourField(const DimensionedField<Type, volMesh>& iField)"
    );

    return Field<Type>::null();
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::patchNeighbourField() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "patchNeighbourField()"
    );

    return Field<Type>::null();
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "valueInternalCoeffs(const tmp<scalarField>&)"
    );

    return Field<Type>::null();
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "valueBoundaryCoeffs(const tmp<scalarField>&)"
    );

    return Field<Type>::null();
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::gradientInternalCoeffs() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "gradientInternalCoeffs()"
    );

    return Field<Type>::null();
}


template<class Type>
tmp<Field<Type> > slicedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    notImplemented
    (
        "slicedFvPatchField<Type>::"
        "gradientBoundaryCoeffs()"
    );

    return Field<Type>::null();
}


template<class Type>
void slicedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
