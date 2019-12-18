/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "CONSTRUCT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::CLASS::t() const
{
    return this->db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(p, iF),
    scalarData_(0.0),
    data_(Zero),
    fieldData_(p.size(), Zero),
    timeVsData_(),
    wordData_("wordDefault"),
    labelData_(-1),
    boolData_(false)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const dictionary& dict
)
:
    PARENT(p, iF),
    scalarData_(dict.lookup<scalar>("scalarData")),
    data_(pTraits<TYPE>(dict.lookup("data"))),
    fieldData_("fieldData", dict, p.size()),
    timeVsData_(Function1<TYPE>::New("timeVsData", dict)),
    wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
    labelData_(-1),
    boolData_(false)
{
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    this->refValue() = FIELD("fieldData", dict, p.size());
    FVPATCHF::operator=(this->refValue());

    PARENT::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    FVPATCHF::operator=
    (
        FIELD("value", dict, p.size())
    );
    this->refValue() = *this;
    */
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    PARENT(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(mapper(ptf.fieldData_)),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf
)
:
    PARENT(ptf),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(ptf, iF),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_, false),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::CLASS::autoMap
(
    const fvPatchFieldMapper& m
)
{
    PARENT::autoMap(m);
    m(fieldData_, fieldData_);
}


template<class Type>
void Foam::CLASS::rmap
(
    const FVPATCHF& ptf,
    const labelList& addr
)
{
    PARENT::rmap(ptf, addr);

    const CLASS& tiptf =
        refCast<const CLASS>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


template<class Type>
void Foam::CLASS::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    PARENT::operator==
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );

    const scalarField& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            "phi"
        );
    this->valueFraction() = 1.0 - pos0(phip);

    PARENT::updateCoeffs();
}


template<class Type>
void Foam::CLASS::write
(
    Ostream& os
) const
{
    FVPATCHF::write(os);
    writeEntry(os, "scalarData", scalarData_);
    writeEntry(os, "data", data_);
    writeEntry(os, "fieldData", fieldData_);
    writeEntry(os, timeVsData_());
    writeEntry(os, "wordData", wordData_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        FVPATCHF,
        CLASS
    );
}

// ************************************************************************* //
