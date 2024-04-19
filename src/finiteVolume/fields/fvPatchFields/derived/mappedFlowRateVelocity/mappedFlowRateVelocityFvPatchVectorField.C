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

#include "mappedFlowRateVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "mappedFvPatchBaseBase.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFlowRateVelocityFvPatchVectorField::
mappedFlowRateVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nbrPhiName_(dict.lookupOrDefault<word>("nbrPhi", "phi")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    mappedPatchBaseBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBaseBase::from::differentPatch
    );
}


Foam::mappedFlowRateVelocityFvPatchVectorField::
mappedFlowRateVelocityFvPatchVectorField
(
    const mappedFlowRateVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::mappedFlowRateVelocityFvPatchVectorField::
mappedFlowRateVelocityFvPatchVectorField
(
    const mappedFlowRateVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFlowRateVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the mapper and the neighbouring patch
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& nbrPatch = mapper.nbrFvPatch();

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField phip
    (
        mapper.fromNeighbour
        (
            nbrPatch.lookupPatchField<surfaceScalarField, scalar>(nbrPhiName_)
        )
    );

    const scalarField Unp(-phip/patch().magSf());

    const vectorField np(patch().nf());

    if (phi.dimensions() == dimVolumetricFlux)
    {
        // volumetric flow-rate
        operator==(np*Unp);
    }
    else if (phi.dimensions() == dimMassFlux)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // mass flow-rate
        operator==(np*Unp/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::mappedFlowRateVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry(os, "nbrPhi", nbrPhiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       mappedFlowRateVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
