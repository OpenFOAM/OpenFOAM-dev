/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "drag.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToConstantDensityFluid.H"
#include "coupledToFluid.H"
#include "LagrangianmSp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(drag, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::drag::addUSup
(
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    const LagrangianSubScalarField& D = this->D(U.mesh());

    const LagrangianSubVectorField& Uc =
        cloud<clouds::carried>().Uc(U.mesh());

    if (eqn.isPsi(U))
    {
        eqn.Su += D*Uc;
        eqn -= Lagrangianm::Sp(D, U);
    }
    else
    {
        eqn += Lagrangianm::Sp(D, Uc);
        eqn.Su -= D*U;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::drag::drag
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    D
    (
        cloud().derivedField<scalar>
        (
            "D",
            *this,
            &drag::calcD
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::drag::addSupFields() const
{
    return wordList({cloud().U.name()});
}


bool Foam::Lagrangian::drag::addsSupToField
(
    const word& fieldName,
    const word& eqnFieldName
) const
{
    return
        fieldName == cloud().U.name()
     && (
            eqnFieldName == cloud().U.name()
         || eqnFieldName == cloud<clouds::coupled>().Uc.name()
        );
}


void Foam::Lagrangian::drag::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToConstantDensityFluid>();

    addUSup(U, eqn);
}


void Foam::Lagrangian::drag::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud
    <
        clouds::coupledToConstantDensityFluid,
        clouds::coupledToFluid
    >();

    addUSup(U, eqn);
}


// ************************************************************************* //
