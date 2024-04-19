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

#include "surfaceNormalFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    refValue_("refValue", iF.dimensions(), dict, p.size())
{
    fvPatchVectorField::operator=(refValue_*patch().nf());
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(mapper(ptf.refValue_))
{
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
    fvPatchVectorField::operator=
    (
        mapper(ptf.refValue_*ptf.patch().nf())
    );
}


Foam::surfaceNormalFixedValueFvPatchVectorField::
surfaceNormalFixedValueFvPatchVectorField
(
    const surfaceNormalFixedValueFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceNormalFixedValueFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const surfaceNormalFixedValueFvPatchVectorField& tiptf =
        refCast<const surfaceNormalFixedValueFvPatchVectorField>(ptf);

    mapper(refValue_, tiptf.refValue_);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const surfaceNormalFixedValueFvPatchVectorField& tiptf =
        refCast<const surfaceNormalFixedValueFvPatchVectorField>(ptf);

    refValue_.reset(tiptf.refValue_);
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchVectorField::operator=(refValue_*patch().nf());
    fvPatchVectorField::updateCoeffs();
}


void Foam::surfaceNormalFixedValueFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "refValue", refValue_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        surfaceNormalFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
