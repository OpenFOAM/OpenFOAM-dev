/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "rotatingPressureInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
calcTangentialVelocity()
{
    const scalar t = this->db().time().timeOutputValue();
    vector om = omega_->value(t);

    vector axisHat = om/mag(om);
    const vectorField tangentialVelocity
    (
        (-om) ^ (patch().Cf() - axisHat*(axisHat & patch().Cf()))
    );

    const vectorField n(patch().nf());
    refValue() = tangentialVelocity - n*(n & tangentialVelocity);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletOutletVelocityFvPatchVectorField(p, iF),
    omega_()
{}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    pressureInletOutletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_().clone().ptr())
{
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    pressureInletOutletVelocityFvPatchVectorField(p, iF, dict),
    omega_(DataEntry<vector>::New("omega", dict))
{
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& rppvf
)
:
    pressureInletOutletVelocityFvPatchVectorField(rppvf),
    omega_(rppvf.omega_().clone().ptr())
{
    calcTangentialVelocity();
}


Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::
rotatingPressureInletOutletVelocityFvPatchVectorField
(
    const rotatingPressureInletOutletVelocityFvPatchVectorField& rppvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletOutletVelocityFvPatchVectorField(rppvf, iF),
    omega_(rppvf.omega_().clone().ptr())
{
    calcTangentialVelocity();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatingPressureInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("phi") << phiName() << token::END_STATEMENT << nl;
    omega_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rotatingPressureInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
