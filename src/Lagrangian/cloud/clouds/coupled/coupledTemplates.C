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

#include "coupled.H"
#include "CarrierField.H"
#include "CarrierEqn.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
const Foam::CarrierField<Type>& Foam::clouds::coupled::carrierField
(
    const VolField<Type>& psi
) const
{
    return carrierField<Type, VolField<Type>>(psi);
}


template<class Type, class ... Args>
const Foam::CarrierField<Type>& Foam::clouds::coupled::carrierField
(
    const Args& ... args
) const
{
    const CarrierField<Type>* lookupPtr =
        carrierFields<Type>().lookupPtr(CarrierField<Type>(args ...).name());

    if (lookupPtr) return *lookupPtr;

    CarrierField<Type>* ptr = new CarrierField<Type>(args ...);
    carrierFields<Type>().insert(ptr->name(), ptr);
    return *ptr;
}


template<class Type>
Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const VolField<Type>& psi
)
{
    CarrierEqn<Type>* lookupPtr =
        carrierEqns<Type>().lookupPtr(psi.name());

    if (lookupPtr) return *lookupPtr;

    CarrierEqn<Type>* ptr = new CarrierEqn<Type>(psi);
    carrierEqns<Type>().insert(psi.name(), ptr);
    return *ptr;
}


template<class Type>
const Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const VolField<Type>& psi
) const
{
    return carrierEqns<Type>()[psi.name()];
}


template<class Type>
Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const CarrierField<Type>& psic
)
{
    return carrierEqn(psic.psi());
}


template<class Type>
const Foam::CarrierEqn<Type>& Foam::clouds::coupled::carrierEqn
(
    const CarrierField<Type>& psic
) const
{
    return carrierEqn(psic.psi());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
