/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "waveAlpha_DimensionedFieldFunction.H"
#include "fvMesh.H"
#include "waveSuperposition.H"
#include "levelSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace DimensionedFieldFunctions
{
    defineTypeNameAndDebug
    (
        waveAlpha,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<scalar, fvMesh>>
        scalarFvMeshDimensionedFieldFunction;
    addToRunTimeSelectionTable
    (
        scalarFvMeshDimensionedFieldFunction,
        waveAlpha,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::waveAlpha::waveAlpha
(
    const dictionary& dict,
    DimensionedField<scalar, fvMesh>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvMesh>>(dict, field),
    liquid_(dict.lookupOrDefault<Switch>("liquid", true))
{}


Foam::DimensionedFieldFunctions::waveAlpha::waveAlpha
(
    const waveAlpha& dff,
    DimensionedField<scalar, fvMesh>& field
)
:
    DimensionedFieldFunction<DimensionedField<scalar, fvMesh>>(dff, field),
    liquid_(dff.liquid_)
{}


Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::scalar, Foam::fvMesh>
    >
>
Foam::DimensionedFieldFunctions::waveAlpha::clone
(
    DimensionedField<scalar, fvMesh>& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedField<scalar, fvMesh>>>
    (
        new waveAlpha(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::waveAlpha::evaluate()
{
    const fvMesh& mesh(this->field_.mesh());
    const scalar t = mesh.time().value();

    const waveSuperposition& waves = waveSuperposition::New(mesh);

    this->field_.primitiveFieldRef() = levelSetFraction
    (
        mesh,
        waves.height(t, mesh.cellCentres()),
        waves.height(t, mesh.points()),
        !liquid_
    );
}


void Foam::DimensionedFieldFunctions::waveAlpha::write
(
    Ostream& os
) const
{
    writeEntryIfDifferent<Switch>(os, "liquid", true, liquid_);
}


// ************************************************************************* //
