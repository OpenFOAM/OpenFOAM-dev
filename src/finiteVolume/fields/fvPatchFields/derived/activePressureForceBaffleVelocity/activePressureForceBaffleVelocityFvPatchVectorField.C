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

#include "activePressureForceBaffleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activePressureForceBaffleVelocityFvPatchVectorField::
activePressureForceBaffleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    orientation_(dict.lookup<label>("orientation")),
    initWallSf_(0),
    initCyclicSf_(0),
    nbrCyclicSf_(0),
    openFraction_(dict.lookup<scalar>("openFraction")),
    openingTime_(dict.lookup<scalar>("openingTime")),
    maxOpenFractionDelta_(dict.lookup<scalar>("maxOpenFractionDelta")),
    curTimeIndex_(-1),
    minThresholdValue_(dict.lookup<scalar>("minThresholdValue")),
    fBased_(readBool(dict.lookup("forceBased"))),
    baffleActivated_(0)
{
    fvPatchVectorField::operator=(Zero);

    if (p.size() > 0)
    {
        initWallSf_ = p.Sf();
        initCyclicSf_ = p.boundaryMesh()[cyclicPatchLabel_].Sf();
        nbrCyclicSf_ =  refCast<const cyclicFvPatch>
        (
            p.boundaryMesh()[cyclicPatchLabel_]
        ).neighbFvPatch().Sf();
    }

    if (dict.found("p"))
    {
        dict.lookup("p") >> pName_;
    }
}


Foam::activePressureForceBaffleVelocityFvPatchVectorField::
activePressureForceBaffleVelocityFvPatchVectorField
(
    const activePressureForceBaffleVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    maxOpenFractionDelta_(ptf.maxOpenFractionDelta_),
    curTimeIndex_(-1),
    minThresholdValue_(ptf.minThresholdValue_),
    fBased_(ptf.fBased_),
    baffleActivated_(ptf.baffleActivated_)
{}


Foam::activePressureForceBaffleVelocityFvPatchVectorField::
activePressureForceBaffleVelocityFvPatchVectorField
(
    const activePressureForceBaffleVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    maxOpenFractionDelta_(ptf.maxOpenFractionDelta_),
    curTimeIndex_(-1),
    minThresholdValue_(ptf.minThresholdValue_),
    fBased_(ptf.fBased_),
    baffleActivated_(ptf.baffleActivated_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activePressureForceBaffleVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    //- Note: cannot map field from cyclic patch anyway so just recalculate
    //  Areas should be consistent when doing map except in case of topo
    //  changes.
    //- Note: we don't want to use Sf here since triggers rebuilding of
    //  fvMesh::S() which will give problems when mapped (since already
    //  on new mesh)

    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::activePressureForceBaffleVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    // See rmap.
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::activePressureForceBaffleVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const volScalarField& p = db().lookupObject<volScalarField>
        (
            pName_
        );

        const fvPatch& cyclicPatch = patch().boundaryMesh()[cyclicPatchLabel_];
        const labelList& cyclicFaceCells = cyclicPatch.patch().faceCells();
        const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();

        const labelList& nbrFaceCells = nbrPatch.patch().faceCells();

        scalar valueDiff = 0;

        if (fBased_)
        {
             // Add this side
            forAll(cyclicFaceCells, facei)
            {
                valueDiff +=p[cyclicFaceCells[facei]]*mag(initCyclicSf_[facei]);
            }

            // Remove other side
            forAll(nbrFaceCells, facei)
            {
                valueDiff -=p[nbrFaceCells[facei]]*mag(initCyclicSf_[facei]);
            }

            Info<< "Force difference = " << valueDiff << endl;
        }
        else // pressure based
        {
            forAll(cyclicFaceCells, facei)
            {
                valueDiff += p[cyclicFaceCells[facei]];
            }

            forAll(nbrFaceCells, facei)
            {
                valueDiff -= p[nbrFaceCells[facei]];
            }

            Info<< "Pressure difference = " << valueDiff << endl;
        }

        if ((mag(valueDiff) > mag(minThresholdValue_)) || baffleActivated_)
        {
            openFraction_ =
                max(
                    min(
                        openFraction_
                      + min
                        (
                          this->db().time().deltaT().value()/openingTime_,
                          maxOpenFractionDelta_
                        )*(orientation_),
                        1 - 1e-6
                        ),
                    1e-6
                    );

             baffleActivated_ = true;
        }
        else
        {
            openFraction_ = max(min(1 - 1e-6, openFraction_), 1e-6);
        }

        Info<< "Open fraction = " << openFraction_ << endl;

        vectorField::subField Sfw = patch().patch().faceAreas();
        vectorField newSfw((1 - openFraction_)*initWallSf_);
        forAll(Sfw, facei)
        {
            Sfw[facei] = newSfw[facei];
        }
        const_cast<scalarField&>(patch().magSf()) = mag(patch().Sf());

        // Update owner side of cyclic
        const_cast<vectorField&>(cyclicPatch.Sf()) =
            openFraction_*initCyclicSf_;

        const_cast<scalarField&>(cyclicPatch.magSf()) =
            mag(cyclicPatch.Sf());

        // Update neighbour side of cyclic
        const_cast<vectorField&>(nbrPatch.Sf()) =
            openFraction_*nbrCyclicSf_;

        const_cast<scalarField&>(nbrPatch.magSf()) =
            mag(nbrPatch.Sf());

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::activePressureForceBaffleVelocityFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry(os, "cyclicPatch", cyclicPatchName_);
    writeEntry(os, "orientation", orientation_);
    writeEntry(os, "openingTime", openingTime_);
    writeEntry(os, "maxOpenFractionDelta", maxOpenFractionDelta_);
    writeEntry(os, "openFraction", openFraction_);
    writeEntry(os, "minThresholdValue", minThresholdValue_);
    writeEntry(os, "forceBased", fBased_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        activePressureForceBaffleVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
