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

#include "FluidMulticomponentThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::FluidMulticomponentThermo<BaseThermo>::FluidMulticomponentThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::FluidMulticomponentThermo<BaseThermo>::~FluidMulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::scalar Foam::FluidMulticomponentThermo<BaseThermo>::mui
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).mu(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::FluidMulticomponentThermo<BaseThermo>::mui
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return this->volScalarFieldPropertyi
    (
        "mu",
        dimMass/dimLength/dimTime,
        &BaseThermo::mixtureType::thermoType::mu,
        speciei,
        p,
        T
    );
}


// ************************************************************************* //
