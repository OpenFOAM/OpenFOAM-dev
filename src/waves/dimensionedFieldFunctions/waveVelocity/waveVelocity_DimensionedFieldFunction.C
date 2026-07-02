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

#include "waveVelocity_DimensionedFieldFunction.H"
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
        waveVelocity,
        0
    );

    typedef DimensionedFieldFunction<DimensionedField<vector, fvMesh>>
        DimensionedFieldFunctionVectorFvMesh;
    addToRunTimeSelectionTable
    (
        DimensionedFieldFunctionVectorFvMesh,
        waveVelocity,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DimensionedFieldFunctions::waveVelocity::waveVelocity
(
    const dictionary& dict,
    DimensionedField<vector, fvMesh>& field
)
:
    DimensionedFieldFunction<DimensionedField<vector, fvMesh>>(dict, field)
{}


Foam::DimensionedFieldFunctions::waveVelocity::waveVelocity
(
    const waveVelocity& dff,
    DimensionedField<vector, fvMesh>& field
)
:
    DimensionedFieldFunction<DimensionedField<vector, fvMesh>>(dff, field)
{}


Foam::autoPtr
<
    Foam::DimensionedFieldFunction
    <
        Foam::DimensionedField<Foam::vector, Foam::fvMesh>
    >
>
Foam::DimensionedFieldFunctions::waveVelocity::clone
(
    DimensionedField<vector, fvMesh>& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedField<vector, fvMesh>>>
    (
        new waveVelocity(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DimensionedFieldFunctions::waveVelocity::evaluate()
{
    const fvMesh& mesh(this->field_.mesh());
    const scalar t = mesh.time().value();

    const waveSuperposition& waves = waveSuperposition::New(mesh);

    const pointField& ccs = mesh.cellCentres();
    const pointField& pts = mesh.points();

    this->field_.primitiveFieldRef() = levelSetAverage
    (
        mesh,
        waves.height(t, ccs),
        waves.height(t, pts),
        waves.UGas(t, ccs)(),
        waves.UGas(t, pts)(),
        waves.ULiquid(t, ccs)(),
        waves.ULiquid(t, pts)()
    );
}


void Foam::DimensionedFieldFunctions::waveVelocity::write
(
    Ostream& os
) const
{}


// ************************************************************************* //
