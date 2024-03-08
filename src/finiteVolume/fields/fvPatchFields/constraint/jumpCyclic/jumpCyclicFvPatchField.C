/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "jumpCyclicFvPatchField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * /

template<class Type>
template<class ... Cmpts>
void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrixCmpts
(
    Field<Type>& result,
    const Field<Type>& psi,
    const scalarField& coeffs,
    const Pstream::commsTypes,
    const Cmpts ... cmpts
) const
{
    const labelUList& faceCells = this->patch().faceCells();
    const labelUList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    Field<Type> nbrPf(psi, nbrFaceCells);

    // Only apply the jump to the original field
    if (&psi == &this->primitiveField())
    {
        nbrPf += jump();
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(nbrPf, cmpts ...);

    // Multiply the field by coefficients and add into the result
    forAll(faceCells, facei)
    {
        result[faceCells[facei]] -= coeffs[facei]*nbrPf[facei];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const jumpCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::jumpCyclicFvPatchField<Type>::jumpCyclicFvPatchField
(
    const jumpCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::jumpCyclicFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    const VolField<Type>& vf =
        static_cast<const VolField<Type>&>(this->internalField());

    const labelUList& nbrFaceCells =
        this->cyclicPatch().neighbFvPatch().faceCells();

    return this->transform().transform(Field<Type>(vf, nbrFaceCells)) + jump();
}


template<class Type>
void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psi,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    NotImplemented;
}


template<class Type>
void Foam::jumpCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psi,
    const scalarField& coeffs,
    const Pstream::commsTypes comms
) const
{
    updateInterfaceMatrixCmpts(result, psi, coeffs, comms);
}


// ************************************************************************* //
