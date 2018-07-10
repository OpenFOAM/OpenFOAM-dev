/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
#include "wavePressureFvPatchScalarField.H"
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
    directionMixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    pName_("p"),
    inletOutlet_(true),
    waves_(db()),
    faceCellSubset_(nullptr),
    faceCellSubsetTimeIndex_(-1)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    inletOutlet_(dict.lookupOrDefault<Switch>("inletOutlet", true)),
    waves_(db(), dict),
    faceCellSubset_(nullptr),
    faceCellSubsetTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    inletOutlet_(ptf.inletOutlet_),
    waves_(ptf.waves_),
    faceCellSubset_(nullptr),
    faceCellSubsetTimeIndex_(-1)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    inletOutlet_(ptf.inletOutlet_),
    waves_(ptf.waves_),
    faceCellSubset_(nullptr),
    faceCellSubsetTimeIndex_(-1)
{}


Foam::waveVelocityFvPatchVectorField::waveVelocityFvPatchVectorField
(
    const waveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    inletOutlet_(ptf.inletOutlet_),
    waves_(ptf.waves_),
    faceCellSubset_(nullptr),
    faceCellSubsetTimeIndex_(-1)
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


Foam::tmp<Foam::vectorField> Foam::waveVelocityFvPatchVectorField::U() const
{
    const scalar t = db().time().timeOutputValue();

    return
        levelSetAverage
        (
            patch(),
            waves_.height(t, patch().Cf()),
            waves_.height(t, patch().patch().localPoints()),
            waves_.UGas(t, patch().Cf())(),
            waves_.UGas(t, patch().patch().localPoints())(),
            waves_.ULiquid(t, patch().Cf())(),
            waves_.ULiquid(t, patch().patch().localPoints())()
        );
}


Foam::tmp<Foam::vectorField> Foam::waveVelocityFvPatchVectorField::Un() const
{
    const scalar t = db().time().timeOutputValue();

    const fvMeshSubset& subset = faceCellSubset();
    const fvMesh& meshs = subset.subMesh();
    const label patchis = findIndex(subset.patchMap(), patch().index());

    const vectorField Us
    (
        levelSetAverage
        (
            meshs,
            waves_.height(t, meshs.cellCentres())(),
            waves_.height(t, meshs.points())(),
            waves_.UGas(t, meshs.cellCentres())(),
            waves_.UGas(t, meshs.points())(),
            waves_.ULiquid(t, meshs.cellCentres())(),
            waves_.ULiquid(t, meshs.points())()
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

    const fvPatchScalarField& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    if (isA<wavePressureFvPatchScalarField>(pp))
    {
        const vectorField U(this->U()), Un(this->Un());
        const scalarField out(pos0(U & patch().Sf()));

        // Where inflow, set all velocity components to values specified by the
        // wave model. Where outflow, set the tangential values and the normal
        // gradient.
        valueFraction() = symmTensor::I - out*sqr(patch().nf());
        refValue() = U;
        refGrad() = (U - Un)*patch().deltaCoeffs();
    }
    else
    {
        const vectorField U(this->U());

        if (inletOutlet_)
        {
            const scalarField& phip =
                patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
            const scalarField out(pos0(phip));

            // Where inflow, fix all velocity components to values specified by
            // the wave model.
            refValue() = (1 - out)*U;
            valueFraction() = (1 - out)*symmTensor::I;

            // Where outflow, set the normal component of the velocity to a
            // value consistent with phi, but scale it to get the volumetric
            // flow rate specified by the wave model. Tangential components are
            // extrapolated.
            const scalar QPhip = gSum(out*phip);
            const scalar QWave = gSum(out*(U & patch().Sf()));
            const vectorField nBySf(patch().Sf()/sqr(patch().magSf()));
            if (QPhip > vSmall)
            {
                refValue() += out*(QWave/QPhip)*phip*nBySf;
            }
            else
            {
                refValue() += out*QWave*nBySf;
            }
            valueFraction() += out*sqr(patch().nf());
        }
        else
        {
            refValue() = U;
            valueFraction() = symmTensor::I;
        }
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::waveVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<Switch>(os, "inletOutlet", true, inletOutlet_);
    waves_.write(os);
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
