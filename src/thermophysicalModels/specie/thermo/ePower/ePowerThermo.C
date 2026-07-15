/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "ePowerThermo.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::ePowerThermo<EquationOfState>::ePowerThermo
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    EquationOfState(name, dict),
    c0_(subDict.lookup<scalar>("c0", dimensions::specificHeatCapacity)),
    n0_(subDict.lookup<scalar>("n0", dimless)),
    Tref_(subDict.lookup<scalar>("Tref", dimensions::temperature)),
    hf_
    (
        subDict.lookupBackwardsCompatible<scalar>
        (
            {"hf", "Hf"},
            dimensions::specificEnergy
        )
    )
{}


template<class EquationOfState>
Foam::ePowerThermo<EquationOfState>::ePowerThermo
(
    const word& name,
    const dictionary& dict
)
:
    ePowerThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::ePowerThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    writeEntry
    (
        os,
        "thermodynamics",
        dictionary::entries("c0", c0_, "n0", n0_, "Tref", Tref_, "hf", hf_)
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ePowerThermo<EquationOfState>& et
)
{
    et.write(os);
    return os;
}


// ************************************************************************* //
