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

#include "mappedVelocityFluxFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedVelocityFluxFvPatchField::mappedVelocityFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi")
{}


Foam::mappedVelocityFluxFvPatchField::mappedVelocityFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    mappedPatchBase::validateMapForField
    (
        *this,
        dict,
        mappedPatchBase::from::differentPatch
    );
}


Foam::mappedVelocityFluxFvPatchField::mappedVelocityFluxFvPatchField
(
    const mappedVelocityFluxFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


Foam::mappedVelocityFluxFvPatchField::mappedVelocityFluxFvPatchField
(
    const mappedVelocityFluxFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedVelocityFluxFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    const mappedPatchBase& mapper = mappedPatchBase::getMap(patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper.nbrMesh());
    const label nbrPatchi = mapper.nbrPolyPatch().index();

    const volVectorField& UField =
        refCast<const volVectorField>(internalField());
    surfaceScalarField& phiField =
        nbrMesh.lookupObjectRef<surfaceScalarField>(phiName_);

    operator==(mapper.distribute(UField.boundaryField()[nbrPatchi]));
    phiField.boundaryFieldRef()[patch().index()] ==
        mapper.distribute(phiField.boundaryField()[nbrPatchi]);

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::mappedVelocityFluxFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mappedVelocityFluxFvPatchField
    );
}


// ************************************************************************* //
