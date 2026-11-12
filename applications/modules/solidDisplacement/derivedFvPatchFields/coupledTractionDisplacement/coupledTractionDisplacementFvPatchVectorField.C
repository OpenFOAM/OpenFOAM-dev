/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "coupledTractionDisplacementFvPatchVectorField.H"
#include "solidDisplacementThermo.H"
#include "mappedFvPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledTractionDisplacementFvPatchVectorField::
coupledTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    rhoFluid_(dict.lookupOrDefault<scalar>("rhoFluid", dimensions::density, 1))
{}


Foam::coupledTractionDisplacementFvPatchVectorField::
coupledTractionDisplacementFvPatchVectorField
(
    const coupledTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper),
    rhoFluid_(tdpvf.rhoFluid_)
{}


Foam::coupledTractionDisplacementFvPatchVectorField::
coupledTractionDisplacementFvPatchVectorField
(
    const coupledTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, iF),
    rhoFluid_(tdpvf.rhoFluid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledTractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the mapper and the neighbouring patch
    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& patchNbr = mapper.nbrFvPatch();

    const scalarField pressure
    (
        mapper.fromNeighbour
        (
            patchNbr.lookupPatchField<volScalarField, scalar>("p")
        )
    );

    tractionDisplacementFvPatchVectorField::updateCoeffs
    (
        rhoFluid_*pressure
    );

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::coupledTractionDisplacementFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "rhoFluid", rhoFluid_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        coupledTractionDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
