/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "symmetryPlaneFvPatchField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(p, iF),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p))
{}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchField<Type>(ptf, p, iF, mapper),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p))
{
    if (!isType<symmetryPlaneFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchField<Type>(p, iF, dict),
    symmetryPlanePatch_(refCast<const symmetryPlaneFvPatch>(p))
{
    if (!isType<symmetryPlaneFvPatch>(p))
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
}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf
)
:
    basicSymmetryFvPatchField<Type>(ptf),
    symmetryPlanePatch_(ptf.symmetryPlanePatch_)
{}


template<class Type>
Foam::symmetryPlaneFvPatchField<Type>::symmetryPlaneFvPatchField
(
    const symmetryPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(ptf, iF),
    symmetryPlanePatch_(ptf.symmetryPlanePatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::symmetryPlaneFvPatchField<Type>::snGrad() const
{
    vector nHat(symmetryPlanePatch_.n());

    const Field<Type> iF(this->patchInternalField());

    return
        (transform(I - 2.0*sqr(nHat), iF) - iF)
       *(this->patch().deltaCoeffs()/2.0);
}


template<class Type>
void Foam::symmetryPlaneFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vector nHat(symmetryPlanePatch_.n());

    const Field<Type> iF(this->patchInternalField());

    Field<Type>::operator=
    (
        (iF + transform(I - 2.0*sqr(nHat), iF))/2.0
    );

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::symmetryPlaneFvPatchField<Type>::snGradTransformDiag() const
{
    vector nHat(symmetryPlanePatch_.n());

    const vector diag
    (
        mag(nHat.component(vector::X)),
        mag(nHat.component(vector::Y)),
        mag(nHat.component(vector::Z))
    );

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            this->size(),
            transformMask<Type>
            (
                // pow<vector, pTraits<Type>::rank>(diag)
                pow
                (
                    diag,
                    pTraits<typename powProduct<vector, pTraits<Type>::rank>
                    ::type>::zero
                )
            )
        )
    );
}


// ************************************************************************* //
