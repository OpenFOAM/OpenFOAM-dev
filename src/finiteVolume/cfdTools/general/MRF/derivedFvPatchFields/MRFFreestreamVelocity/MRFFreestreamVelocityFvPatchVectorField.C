/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "MRFFreestreamVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFFreestreamVelocityFvPatchVectorField::
MRFFreestreamVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    freestreamVelocityFvPatchVectorField(p, iF, dict),
    MRFPatchField(dict),
    freestreamValue0_(dict.lookup("freestreamValue0"))
{}


Foam::MRFFreestreamVelocityFvPatchVectorField::
MRFFreestreamVelocityFvPatchVectorField
(
    const MRFFreestreamVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    freestreamVelocityFvPatchVectorField(pvf, p, iF, mapper),
    MRFPatchField(pvf),
    freestreamValue0_(pvf.freestreamValue0_)
{}


Foam::MRFFreestreamVelocityFvPatchVectorField::
MRFFreestreamVelocityFvPatchVectorField
(
    const MRFFreestreamVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    freestreamVelocityFvPatchVectorField(pvf, iF),
    MRFPatchField(pvf),
    freestreamValue0_(pvf.freestreamValue0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFFreestreamVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar time = this->db().time().value();
    const vector Omega = MRFzone(db()).Omega();
    const scalar omega = mag(Omega);
    const vector axis(Omega/omega);
    const scalar theta = time*omega;

    refValue() =
        cos(theta)*freestreamValue0_
      - sin(theta)*(axis ^ freestreamValue0_);

    freestreamVelocityFvPatchVectorField::updateCoeffs();
}


void Foam::MRFFreestreamVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "freestreamValue0", freestreamValue0_);
    writeEntry(os, "freestreamValue", freestreamValue());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        MRFFreestreamVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
