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

#include "thermal.H"
#include "CloudTypes.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class ... Thermos>
bool Foam::clouds::thermal::isThermo() const
{
    return CloudTypes<Thermo, Thermos ...>::isA(autoPtr::operator()());
}


template<class Thermo, class ... Thermos>
void Foam::clouds::thermal::assertThermo() const
{
    if (!isThermo<Thermo, Thermos ...>())
    {
        FatalErrorInFunction
            << "The cloud '"
            << cloud_.mesh().name() << "' requires a thermodynamic model "
            << "derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '" << autoPtr::operator()().type()
            << "'" << exit(FatalError);
    }
}


template<class Thermo, class ... Thermos>
void Foam::clouds::thermal::assertThermo(const LagrangianModel& model) const
{
    if (!isThermo<Thermo, Thermos ...>())
    {
        FatalErrorInFunction
            << "The Larangian model '" << model.name() << "' of cloud '"
            << cloud_.mesh().name() << "' requires a thermodynamic model "
            << "derived from "
            << CloudTypes<Thermo, Thermos ...>::typesString("or").c_str()
            << ", rather than '" << autoPtr::operator()().type()
            << "'" << exit(FatalError);
    }
}


template<class Thermo, class ... Args>
const Thermo& Foam::clouds::thermal::thermo(const Args& ... args) const
{
    assertThermo<Thermo>(args ...);
    return refCast<const Thermo>(autoPtr::operator()());
}


template<class Thermo, class ... Args>
Thermo& Foam::clouds::thermal::thermo(const Args& ... args)
{
    assertThermo<Thermo>(args ...);
    return refCast<Thermo>(autoPtr::operator()());
}


// ************************************************************************* //
