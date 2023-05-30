/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "heLiquidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class HeThermo>
Foam::heLiquidThermo<HeThermo>::heLiquidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    HeThermo(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class HeThermo>
Foam::heLiquidThermo<HeThermo>::~heLiquidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heLiquidThermo<HeThermo>::sigma() const
{
    return this->volScalarFieldProperty
    (
        "sigma",
        dimForce/dimLength,
        &HeThermo::mixtureType::thermoMixture,
        &HeThermo::mixtureType::thermoMixtureType::sigma,
        this->p_,
        this->T_
    );
}


// ************************************************************************* //
