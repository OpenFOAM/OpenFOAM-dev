/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "polyFacesFvsPatchLabelField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyFacesFvsPatchLabelField::init()
{
    labelField::operator=(identityMap(patch().size()) + patch().start());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(p, iF)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchLabelField(p, iF, dict, false)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const polyFacesFvsPatchLabelField& ptf,
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchLabelField(ptf, p, iF, mapper, false)
{
    init();
}


Foam::polyFacesFvsPatchLabelField::polyFacesFvsPatchLabelField
(
    const polyFacesFvsPatchLabelField& ptf,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(ptf, iF)
{
    init();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFvsPatchTypeField
    (
        fvsPatchLabelField,
        polyFacesFvsPatchLabelField
    );
}


// ************************************************************************* //
