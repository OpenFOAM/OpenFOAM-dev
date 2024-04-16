/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "PhysicalPropertiesThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermoType>
Foam::PhysicalPropertiesThermo<BasicThermoType>::PhysicalPropertiesThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    physicalProperties(mesh, phaseName),
    BasicThermoType(*this, mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermoType>
Foam::PhysicalPropertiesThermo<BasicThermoType>::~PhysicalPropertiesThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasicThermoType>
const Foam::IOdictionary&
Foam::PhysicalPropertiesThermo<BasicThermoType>::properties() const
{
    return *this;
}


template<class BasicThermoType>
Foam::IOdictionary&
Foam::PhysicalPropertiesThermo<BasicThermoType>::properties()
{
    return *this;
}


template<class BasicThermoType>
bool Foam::PhysicalPropertiesThermo<BasicThermoType>::read()
{
    if (physicalProperties::read())
    {
        BasicThermoType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
