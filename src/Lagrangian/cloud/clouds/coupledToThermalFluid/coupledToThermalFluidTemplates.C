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

#include "coupledToThermalFluid.H"
#include "CloudTypes.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class ... Thermos>
bool Foam::clouds::coupledToThermalFluid::isThermoc() const
{
    if (!autoPtr<fluidCarrierThermo>::valid()) return false;

    return
        CloudTypes<Thermo, Thermos ...>::isA
        (
            autoPtr<fluidCarrierThermo>::operator()()
        );
}


template<class Thermo, class ... Thermos>
bool Foam::clouds::coupledToThermalFluid::isThermocs() const
{
    if (!autoPtr<fluidSurfaceThermo>::valid()) return false;

    return
        CloudTypes<Thermo, Thermos ...>::isA
        (
            autoPtr<fluidSurfaceThermo>::operator()()
        );
}


template<class Thermo, class ... Thermos>
bool Foam::clouds::coupledToThermalFluid::isThermocPhase() const
{
    if (!autoPtr<carrierThermo>::valid()) return false;

    return
        CloudTypes<Thermo, Thermos ...>::isA
        (
            autoPtr<carrierThermo>::operator()()
        );
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermoc() const
{
    if (isThermoc<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The cloud '" << cloud_.mesh().name()
        << "' requires the carrier to have a thermodynamic model";

    if (autoPtr<fluidCarrierThermo>::valid())
    {
        FatalError
            << " derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<fluidCarrierThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermocs() const
{
    if (isThermocs<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The cloud '" << cloud_.mesh().name()
        << "' requires the carrier to have a surface thermodynamic model";

    if (autoPtr<fluidSurfaceThermo>::valid())
    {
        FatalError
            << " derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<fluidSurfaceThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermocPhase() const
{
    if (isThermocPhase<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The cloud '" << cloud_.mesh().name() << "' requires the "
        << "corresponding Eulerian phase have a thermodynamic model";

    if (autoPtr<carrierThermo>::valid())
    {
        FatalError
            << " derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<carrierThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermoc
(
    const LagrangianModel& model
) const
{
    if (isThermoc<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The Lagrangian model '" << model.name() << "' of cloud '"
        << cloud_.mesh().name() << "' requires the carrier to have a "
        << "thermodynamic model";

    if (autoPtr<fluidCarrierThermo>::valid())
    {
        FatalError
            << "derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<fluidCarrierThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermocs
(
    const LagrangianModel& model
) const
{
    if (isThermoc<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The Lagrangian model '" << model.name() << "' of cloud '"
        << cloud_.mesh().name() << "' requires the carrier to have a "
        << "surface thermodynamic model";

    if (autoPtr<fluidSurfaceThermo>::valid())
    {
        FatalError
            << "derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<fluidSurfaceThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Thermos>
void Foam::clouds::coupledToThermalFluid::assertThermocPhase
(
    const LagrangianModel& model
) const
{
    if (isThermocPhase<Thermo, Thermos ...>()) return;

    FatalErrorInFunction
        << "The Lagrangian model '" << model.name() << "' of cloud '"
        << cloud_.mesh().name() << "' requires the corresponding Eulerian "
        << "phase to have a thermodynamic model";

    if (autoPtr<carrierThermo>::valid())
    {
        FatalError
            << "derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '"
            << autoPtr<carrierThermo>::operator()().type() << "'";
    }

    FatalError << exit(FatalError);
}


template<class Thermo, class ... Args>
const Thermo& Foam::clouds::coupledToThermalFluid::thermoc
(
    const Args& ... args
) const
{
    assertThermoc<Thermo>(args ...);

    return refCast<const Thermo>(autoPtr<fluidCarrierThermo>::operator()());
}


template<class Thermo, class ... Args>
const Thermo& Foam::clouds::coupledToThermalFluid::thermocs
(
    const Args& ... args
) const
{
    assertThermocs<Thermo>(args ...);

    return refCast<const Thermo>(autoPtr<fluidSurfaceThermo>::operator()());
}


template<class Thermo, class ... Args>
Thermo& Foam::clouds::coupledToThermalFluid::thermocsRef
(
    const Args& ... args
)
{
    assertThermocs<Thermo>(args ...);

    return refCast<Thermo>(autoPtr<fluidSurfaceThermo>::operator()());
}


template<class Thermo, class ... Args>
const Thermo& Foam::clouds::coupledToThermalFluid::thermocPhase
(
    const Args& ... args
) const
{
    assertThermocPhase<Thermo>(args ...);

    return refCast<const Thermo>(autoPtr<carrierThermo>::operator()());
}


// ************************************************************************* //
