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

#include "eTabulatedThermo.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    EquationOfState(name, dict),
    hf_
    (
        subDict.lookupBackwardsCompatible<scalar>
        (
            {"hf", "Hf"},
            dimEnergy/dimMass
        )
    ),
    sf_
    (
        subDict.lookupBackwardsCompatible<scalar>
        (
            {"sf", "Sf"},
            dimEnergy/dimTemperature/dimMass
        )
    ),
    es_
    (
        "es",
        {dimPressure, dimTemperature, dimEnergy/dimMass},
        subDict.subDictBackwardsCompatible({"es", "Es"})
    ),
    Cp_
    (
        "Cp",
        {dimPressure, dimTemperature, dimSpecificHeatCapacity},
        subDict.subDict("Cp")
    ),
    Cv_
    (
        "Cv",
        {dimPressure, dimTemperature, dimSpecificHeatCapacity},
        subDict.subDict("Cv")
    )
{}


template<class EquationOfState>
Foam::eTabulatedThermo<EquationOfState>::eTabulatedThermo
(
    const word& name,
    const dictionary& dict
)
:
    eTabulatedThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eTabulatedThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    const delimitDictionary delimit(os, "thermodynamics");
    writeEntry(os, "hf", hf_);
    writeEntry(os, "sf", sf_);
    {
        const delimitDictionary delimitEs(os, "es");
        es_.write(os, {dimTemperature, dimSpecificHeatCapacity});
    }
    {
        const delimitDictionary delimitCp(os, "Cp");
        Cp_.write(os, {dimTemperature, dimSpecificHeatCapacity});
    }
    {
        const delimitDictionary delimitCv(os, "Cv");
        Cv_.write(os, {dimTemperature, dimSpecificHeatCapacity});
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eTabulatedThermo<EquationOfState>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
