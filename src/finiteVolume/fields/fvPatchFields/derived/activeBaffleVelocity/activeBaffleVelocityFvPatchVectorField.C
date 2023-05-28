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

#include "activeBaffleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
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
    initWallSf_(p.Sf()),
    initCyclicSf_(p.boundaryMesh()[cyclicPatchLabel_].Sf()),
    nbrCyclicSf_
    (
        refCast<const cyclicFvPatch>
        (
            p.boundaryMesh()[cyclicPatchLabel_]
        ).neighbFvPatch().Sf()
    ),
    openFraction_(dict.lookup<scalar>("openFraction")),
    openingTime_(dict.lookup<scalar>("openingTime")),
    maxOpenFractionDelta_(dict.lookup<scalar>("maxOpenFractionDelta")),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(Zero);
}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const activeBaffleVelocityFvPatchVectorField& ptf,
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
    curTimeIndex_(-1)
{}


Foam::activeBaffleVelocityFvPatchVectorField::
activeBaffleVelocityFvPatchVectorField
(
    const activeBaffleVelocityFvPatchVectorField& ptf,
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
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::activeBaffleVelocityFvPatchVectorField::map
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


void Foam::activeBaffleVelocityFvPatchVectorField::reset
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


void Foam::activeBaffleVelocityFvPatchVectorField::updateCoeffs()
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

        scalar forceDiff = 0;

        // Add this side
        forAll(cyclicFaceCells, facei)
        {
            forceDiff += p[cyclicFaceCells[facei]]*mag(initCyclicSf_[facei]);
        }

        // Remove other side
        forAll(nbrFaceCells, facei)
        {
            forceDiff -= p[nbrFaceCells[facei]]*mag(nbrCyclicSf_[facei]);
        }

        openFraction_ =
            max
            (
                min
                (
                    openFraction_
                  + min
                    (
                        this->db().time().deltaTValue()/openingTime_,
                        maxOpenFractionDelta_
                    )
                   *(orientation_*sign(forceDiff)),
                    1 - 1e-6
                ),
                1e-6
            );

        Info<< "openFraction = " << openFraction_ << endl;

        vectorField::subField Sfw = this->patch().patch().faceAreas();
        const vectorField newSfw((1 - openFraction_)*initWallSf_);
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


void Foam::activeBaffleVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry(os, "cyclicPatch", cyclicPatchName_);
    writeEntry(os, "orientation", orientation_);
    writeEntry(os, "openingTime", openingTime_);
    writeEntry(os, "maxOpenFractionDelta", maxOpenFractionDelta_);
    writeEntry(os, "openFraction", openFraction_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        activeBaffleVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
