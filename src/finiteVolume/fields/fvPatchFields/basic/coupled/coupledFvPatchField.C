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

#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF, f)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(p, iF, dict, valueRequired)
{}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(p)),
    fvPatchField<Type>(ptf, p, iF, mapper, false)
{
    mapper(*this, ptf, pTraits<Type>::nan);
}


template<class Type>
Foam::coupledFvPatchField<Type>::coupledFvPatchField
(
    const coupledFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf.patch())),
    fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    mapper(*this, ptf, pTraits<Type>::nan);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coupledFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    return
        deltaCoeffs
       *(this->patchNeighbourField() - this->patchInternalField());
}


template<class Type>
void Foam::coupledFvPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void Foam::coupledFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*this->patchNeighbourField()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -Type(pTraits<Type>::one)*deltaCoeffs;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientInternalCoeffs() const
{
    NotImplemented;
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -this->gradientInternalCoeffs(deltaCoeffs);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    NotImplemented;
    return -this->gradientInternalCoeffs();
}


template<class Type>
void Foam::coupledFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
