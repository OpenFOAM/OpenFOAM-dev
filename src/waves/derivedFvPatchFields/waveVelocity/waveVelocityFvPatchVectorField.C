/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "waveVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "levelSet.H"
#include "volFields.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueInletOutletFvPatchField<vector>(p, iF)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueInletOutletFvPatchField<vector>(p, iF, dict, false)
{
    if (dict.found("value"))
    {
        fixedValueInletOutletFvPatchField<vector>::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fixedValueInletOutletFvPatchField<vector>::operator==
        (
            patchInternalField()
        );
    }
}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueInletOutletFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueInletOutletFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMeshSubset&
Foam::waveVelocityFvPatchVectorField::faceCellSubset() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label timeIndex = mesh.time().timeIndex();

    if
    (
        !faceCellSubset_.valid()
     || (mesh.changing() && faceCellSubsetTimeIndex_ != timeIndex)
    )
    {
        faceCellSubset_.reset(new fvMeshSubset(mesh));
        faceCellSubset_->setCellSubset(patch().faceCells());
        faceCellSubsetTimeIndex_ = timeIndex;

        // Ask for the tetBasePtIs to trigger all processors to build them.
        // Without this, processors that do not contain this patch will
        // generate a comms mismatch.
        faceCellSubset_->subMesh().tetBasePtIs();
    }

    return faceCellSubset_();
}


Foam::tmp<Foam::vectorField>
Foam::waveVelocityFvPatchVectorField::U(const scalar t) const
{
    const waveSuperposition& waves = waveSuperposition::New(db());

    return
        levelSetAverage
        (
            patch(),
            waves.height(t, patch().Cf()),
            waves.height(t, patch().patch().localPoints()),
            waves.UGas(t, patch().Cf())(),
            waves.UGas(t, patch().patch().localPoints())(),
            waves.ULiquid(t, patch().Cf())(),
            waves.ULiquid(t, patch().patch().localPoints())()
        );
}


Foam::tmp<Foam::vectorField>
Foam::waveVelocityFvPatchVectorField::Un(const scalar t) const
{
    const waveSuperposition& waves = waveSuperposition::New(db());

    const fvMeshSubset& subset = faceCellSubset();
    const fvMesh& meshs = subset.subMesh();
    const label patchis = findIndex(subset.patchMap(), patch().index());

    const vectorField Us
    (
        levelSetAverage
        (
            meshs,
            waves.height(t, meshs.cellCentres())(),
            waves.height(t, meshs.points())(),
            waves.UGas(t, meshs.cellCentres())(),
            waves.UGas(t, meshs.points())(),
            waves.ULiquid(t, meshs.cellCentres())(),
            waves.ULiquid(t, meshs.points())()
        )
    );

    tmp<vectorField> tResult(new vectorField(patch().size()));
    vectorField& result = tResult.ref();

    if (patchis != -1)
    {
        forAll(meshs.boundary()[patchis], is)
        {
            const label fs = is + meshs.boundary()[patchis].patch().start();
            const label cs = meshs.boundary()[patchis].faceCells()[is];
            const label f = subset.faceMap()[fs];
            const label i = patch().patch().whichFace(f);
            result[i] = Us[cs];
        }
    }

    return tResult;
}


void Foam::waveVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    operator==(U(t));

    fixedValueInletOutletFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        waveVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
