/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "reboundVelocityLagrangianPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reboundVelocityLagrangianPatchVectorField::
reboundVelocityLagrangianPatchVectorField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    cloudVelocityLagrangianPatchVectorField(p, iIo, dict),
    e_(dict.lookup<scalar>("e", unitless)),
    mu_(dict.lookup<scalar>("mu", unitFraction))
{}


Foam::reboundVelocityLagrangianPatchVectorField::
reboundVelocityLagrangianPatchVectorField
(
    const reboundVelocityLagrangianPatchVectorField& ptf
)
:
    cloudVelocityLagrangianPatchVectorField(ptf),
    e_(ptf.e_),
    mu_(ptf.mu_)
{}


Foam::reboundVelocityLagrangianPatchVectorField::
reboundVelocityLagrangianPatchVectorField
(
    const reboundVelocityLagrangianPatchVectorField& ptf,
    const regIOobject& iIo
)
:
    cloudVelocityLagrangianPatchVectorField(ptf, iIo),
    e_(ptf.e_),
    mu_(ptf.mu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reboundVelocityLagrangianPatchVectorField::evaluate
(
    PstreamBuffers& pBufs,
    const LagrangianScalarInternalDynamicField& fraction
)
{
    cloudVelocityLagrangianPatchVectorField::evaluate(pBufs, fraction);

    const SubField<vector> U = primitiveSubField();

    const vectorField np(patch().mesh().nf<vectorField>(fraction));
    const vectorField Up(patch().mesh().Uf<vectorField>(fraction));

    // !!! Modify Up with tangential velocity from the carrier U boundary
    // condition, to account for moving/not-moving walls, like the lid of a
    // driven cavity

    const scalarField UReln((U - Up) & np);
    const vectorField URelt((U - Up) - np*UReln);

    LagrangianPatchVectorField::operator=(Up - e_*UReln*np + (1 - mu_)*URelt);
}


void Foam::reboundVelocityLagrangianPatchVectorField::write(Ostream& os) const
{
    cloudVelocityLagrangianPatchVectorField::write(os);

    writeEntry(os, "e", e_);
    writeEntry(os, "mu", mu_);
}


Foam::LagrangianState
Foam::reboundVelocityLagrangianPatchVectorField::state() const
{
    return LagrangianState::inCell;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianPatchTypeField
    (
        LagrangianPatchVectorField,
        reboundVelocityLagrangianPatchVectorField
    );
}


// ************************************************************************* //
