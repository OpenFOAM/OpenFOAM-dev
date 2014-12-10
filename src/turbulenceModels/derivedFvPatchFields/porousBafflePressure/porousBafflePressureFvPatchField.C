/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "porousBafflePressureFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::porousBafflePressureFvPatchField<Type>::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    D_(0),
    I_(0),
    length_(0)
{}


template<class Type>
Foam::porousBafflePressureFvPatchField<Type>::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<Type>(ptf, p, iF, mapper),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


template<class Type>
Foam::porousBafflePressureFvPatchField<Type>::porousBafflePressureFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    D_(readScalar(dict.lookup("D"))),
    I_(readScalar(dict.lookup("I"))),
    length_(readScalar(dict.lookup("length")))
{
    fvPatchField<Type>::operator=
    (
        Field<Type>("value", dict, p.size())
    );
}


template<class Type>
Foam::porousBafflePressureFvPatchField<Type>::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<Type>(ptf),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


template<class Type>
Foam::porousBafflePressureFvPatchField<Type>::porousBafflePressureFvPatchField
(
    const porousBafflePressureFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(ptf, iF),
    D_(ptf.D_),
    I_(ptf.I_),
    length_(ptf.length_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::porousBafflePressureFvPatchField<Type>::write(Ostream& os) const
{
    fixedJumpFvPatchField<Type>::write(os);
    os.writeKeyword("D") << D_ << token::END_STATEMENT << nl;
    os.writeKeyword("I") << I_ << token::END_STATEMENT << nl;
    os.writeKeyword("length") << length_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
