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

#include "AttouFerschneider.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(AttouFerschneider, 0);
    addToRunTimeSelectionTable(dragModel, AttouFerschneider, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KGasLiquid
(
    const phaseModel& gas,
    const phaseModel& liquid
) const
{
    const phaseModel& solid = gas.fluid().phases()[solidName_];

    const volScalarField oneMinusGas(max(1 - gas, liquid.residualAlpha()));
    const volScalarField cbrtR
    (
        cbrt(max(solid, solid.residualAlpha())/oneMinusGas)
    );
    const volScalarField magURel(mag(gas.U() - liquid.U()));

    return
        E2_*gas.fluidThermo().mu()*sqr(oneMinusGas/solid.d())*sqr(cbrtR)
       /max(gas, gas.residualAlpha())
      + E2_*gas.rho()*magURel*(1 - gas)/solid.d()*cbrtR;
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KGasSolid
(
    const phaseModel& gas,
    const phaseModel& solid
) const
{
    const volScalarField oneMinusGas(max(1 - gas, solid.residualAlpha()));
    const volScalarField cbrtR
    (
        cbrt(max(solid, solid.residualAlpha())/oneMinusGas)
    );

    return
        E1_*gas.fluidThermo().mu()*sqr(oneMinusGas/solid.d())*sqr(cbrtR)
       /max(gas, gas.residualAlpha())
      + E2_*gas.rho()*mag(gas.U())*(1 - gas)/solid.d()*cbrtR;
}


Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::KLiquidSolid
(
    const phaseModel& liquid,
    const phaseModel& solid
) const
{
    const phaseModel& gas = liquid.fluid().phases()[gasName_];

    return
        E1_*liquid.fluidThermo().mu()
       *sqr(max(solid, solid.residualAlpha())/solid.d())
       /max(liquid, liquid.residualAlpha())
      + E2_*liquid.rho()*mag(gas.U())*solid/solid.d();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::AttouFerschneider::AttouFerschneider
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface),
    gasName_(dict.lookup("gas")),
    liquidName_(dict.lookup("liquid")),
    solidName_(dict.lookup("solid")),
    E1_("E1", dimless, dict),
    E2_("E2", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::AttouFerschneider::~AttouFerschneider()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dragModels::AttouFerschneider::K() const
{
    const phaseModel& gas = interface_.fluid().phases()[gasName_];
    const phaseModel& liquid = interface_.fluid().phases()[liquidName_];
    const phaseModel& solid = interface_.fluid().phases()[solidName_];

    if (interface_.contains(gas) && interface_.contains(liquid))
    {
        return KGasLiquid(gas, liquid);
    }
    if (interface_.contains(gas) && interface_.contains(solid))
    {
        return KGasSolid(gas, solid);
    }
    if (interface_.contains(liquid) && interface_.contains(solid))
    {
        return KLiquidSolid(liquid, solid);
    }

    FatalErrorInFunction
        << "The interface " << interface_.name() << " does not contain two "
        << "out of the gas, liquid and solid phase models."
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::AttouFerschneider::Kf() const
{
    return fvc::interpolate(K());
}


// ************************************************************************* //
