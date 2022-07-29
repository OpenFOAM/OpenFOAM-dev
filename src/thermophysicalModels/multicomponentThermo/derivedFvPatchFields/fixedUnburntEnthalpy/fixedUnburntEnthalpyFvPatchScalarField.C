/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fixedUnburntEnthalpyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "psiuMulticomponentThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedUnburntEnthalpyFvPatchScalarField::
fixedUnburntEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::fixedUnburntEnthalpyFvPatchScalarField::
fixedUnburntEnthalpyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::fixedUnburntEnthalpyFvPatchScalarField::
fixedUnburntEnthalpyFvPatchScalarField
(
    const fixedUnburntEnthalpyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixedUnburntEnthalpyFvPatchScalarField::
fixedUnburntEnthalpyFvPatchScalarField
(
    const fixedUnburntEnthalpyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedUnburntEnthalpyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const psiuMulticomponentThermo& thermo =
        db().lookupObject<psiuMulticomponentThermo>
        (
            physicalProperties::typeName
        );

    const label patchi = patch().index();

    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.Tu().boundaryField()[patchi]);
    Tw.evaluate();
    operator==(thermo.heu(Tw, patchi));

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedUnburntEnthalpyFvPatchScalarField
    );
}

// ************************************************************************* //
