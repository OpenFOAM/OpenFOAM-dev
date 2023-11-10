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

#include "cyclicSlipPointPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicSlipPointPatchField<Type>::cyclicSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(p, iF)
{}


template<class Type>
Foam::cyclicSlipPointPatchField<Type>::cyclicSlipPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    cyclicPointPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::cyclicSlipPointPatchField<Type>::cyclicSlipPointPatchField
(
    const cyclicSlipPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const fieldMapper& mapper
)
:
    cyclicPointPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::cyclicSlipPointPatchField<Type>::cyclicSlipPointPatchField
(
    const cyclicSlipPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    cyclicPointPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::cyclicSlipPointPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    const vectorField& nHat = this->patch().pointNormals();

    tmp<Field<Type>> tvalues =
    (
        (
            this->patchInternalField()
          + transform(I - 2.0*sqr(nHat), this->patchInternalField())
        )/2.0
    );

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInternalField(iF, tvalues());
}


// ************************************************************************* //
