/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "waveAlphaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "levelSet.H"
#include "volFields.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueInletOutletFvPatchScalarField(p, iF, dict, false),
    liquid_(dict.lookupOrDefault<Switch>("liquid", true))
{
    if (dict.found("value"))
    {
        fixedValueInletOutletFvPatchScalarField::operator==
        (
            scalarField("value", iF.dimensions(), dict, p.size())
        );
    }
    else
    {
        fixedValueInletOutletFvPatchScalarField::operator==
        (
            patchInternalField()
        );
    }
}


Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const waveAlphaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueInletOutletFvPatchScalarField(ptf, p, iF, mapper),
    liquid_(ptf.liquid_)
{}


Foam::waveAlphaFvPatchScalarField::waveAlphaFvPatchScalarField
(
    const waveAlphaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueInletOutletFvPatchScalarField(ptf, iF),
    liquid_(ptf.liquid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMeshSubset&
Foam::waveAlphaFvPatchScalarField::faceCellSubset() const
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


Foam::tmp<Foam::scalarField>
Foam::waveAlphaFvPatchScalarField::alpha(const scalar t) const
{
    const waveSuperposition& waves = waveSuperposition::New(db());

    return
        levelSetFraction
        (
            patch(),
            waves.height(t, patch().Cf()),
            waves.height(t, patch().patch().localPoints()),
            !liquid_
        );
}


Foam::tmp<Foam::scalarField>
Foam::waveAlphaFvPatchScalarField::alphan(const scalar t) const
{
    const waveSuperposition& waves = waveSuperposition::New(db());

    const fvMeshSubset& subset = faceCellSubset();
    const fvMesh& meshs = subset.subMesh();
    const label patchis = findIndex(subset.patchMap(), patch().index());

    const scalarField alphas
    (
        levelSetFraction
        (
            meshs,
            waves.height(t, meshs.cellCentres())(),
            waves.height(t, meshs.points())(),
            !liquid_
        )
    );

    tmp<scalarField> tResult(new scalarField(patch().size()));
    scalarField& result = tResult.ref();

    if (patchis != -1)
    {
        forAll(meshs.boundary()[patchis], is)
        {
            const label fs = is + meshs.boundary()[patchis].patch().start();
            const label cs = meshs.boundary()[patchis].faceCells()[is];
            const label f = subset.faceMap()[fs];
            const label i = patch().patch().whichFace(f);
            result[i] = alphas[cs];
        }
    }

    return tResult;
}


void Foam::waveAlphaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(alpha(db().time().value()));

    fixedValueInletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::waveAlphaFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueInletOutletFvPatchScalarField::write(os);
    writeEntryIfDifferent<Switch>(os, "liquid", true, liquid_);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveAlphaFvPatchScalarField
    );
}

// ************************************************************************* //
