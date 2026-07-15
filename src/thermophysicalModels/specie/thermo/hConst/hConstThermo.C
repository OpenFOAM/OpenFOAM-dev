/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "hConstThermo.H"
#include "dictionary.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    EquationOfState(name, dict),
    Cp_(subDict.lookup<scalar>("Cp", dimensions::specificHeatCapacity)),
    hf_
    (
        subDict.lookupBackwardsCompatible<scalar>
        (
            {"hf", "Hf"},
            dimensions::specificEnergy
        )
    ),
    Tref_
    (
        subDict.lookupOrDefault<scalar>
        (
            "Tref",
            dimensions::temperature,
            constant::thermodynamic::Tstd
        )
    ),
    hsRef_
    (
        subDict.lookupOrDefaultBackwardsCompatible<scalar>
        (
            {"hsRef", "Hsref"},
            dimensions::specificEnergy,
            0
        )
    )
{}


template<class EquationOfState>
Foam::hConstThermo<EquationOfState>::hConstThermo
(
    const word& name,
    const dictionary& dict
)
:
    hConstThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hConstThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    const delimitDictionary delimit(os, "thermodynamics");
    writeEntry(os, "Cp", Cp_);
    writeEntry(os, "hf", hf_);
    writeEntryIfDifferent(os, "Tref", constant::thermodynamic::Tstd, Tref_);
    writeEntryIfDifferent(os, "hsRef", scalar(0), hsRef_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hConstThermo<EquationOfState>& ct
)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
