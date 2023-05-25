/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "hydrostaticDisplacementFvPatchVectorField.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hydrostaticDisplacementFvPatchVectorField::
hydrostaticDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    rhoLiquid_(0.0),
    liquidSurfacePressure_(0.0),
    liquidSurfacePoint_(Zero)
{}


Foam::hydrostaticDisplacementFvPatchVectorField::
hydrostaticDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    rhoLiquid_(dict.lookup<scalar>("rhoLiquid")),
    liquidSurfacePressure_(dict.lookup<scalar>("liquidSurfacePressure")),
    liquidSurfacePoint_(dict.lookup("liquidSurfacePoint"))
{}


Foam::hydrostaticDisplacementFvPatchVectorField::
hydrostaticDisplacementFvPatchVectorField
(
    const hydrostaticDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper),
    rhoLiquid_(tdpvf.rhoLiquid_),
    liquidSurfacePressure_(tdpvf.liquidSurfacePressure_),
    liquidSurfacePoint_(tdpvf.liquidSurfacePoint_)
{}


Foam::hydrostaticDisplacementFvPatchVectorField::
hydrostaticDisplacementFvPatchVectorField
(
    const hydrostaticDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, iF),
    rhoLiquid_(tdpvf.rhoLiquid_),
    liquidSurfacePressure_(tdpvf.liquidSurfacePressure_),
    liquidSurfacePoint_(tdpvf.liquidSurfacePoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hydrostaticDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g =
        db().lookupObject<uniformDimensionedVectorField>("g");

    tractionDisplacementFvPatchVectorField::updateCoeffs
    (
        liquidSurfacePressure_
      + rhoLiquid_
       *max(g.value() & (patch().Cf() - liquidSurfacePoint_), scalar(0))
    );
}


void Foam::hydrostaticDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "rhoLiquid", rhoLiquid_);
    writeEntry(os, "liquidSurfacePressure", liquidSurfacePressure_);
    writeEntry(os, "liquidSurfacePoint", liquidSurfacePoint_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        hydrostaticDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
