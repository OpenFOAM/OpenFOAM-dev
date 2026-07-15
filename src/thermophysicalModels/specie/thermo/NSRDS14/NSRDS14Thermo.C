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

#include "NSRDS14Thermo.H"
#include "NSRDS14.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::NSRDS14Thermo<EquationOfState>::NSRDS14Thermo
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
            dimensions::specificEnergy
        )
    ),
    sf_
    (
        subDict.lookupBackwardsCompatible<scalar>
        (
            {"sf", "Sf"},
            dimensions::energy/dimensions::temperature/dimensions::mass
        )
    ),
    Tc_(subDict.lookup<scalar>("Tc", dimensions::temperature)),
    a_(subDict.lookup<scalar>("a", sqrt(dimensions::specificHeatCapacity))),
    b_(subDict.lookup<scalar>("b", dimensions::specificHeatCapacity)),
    c_(subDict.lookup<scalar>("c", sqrt(dimensions::specificHeatCapacity))),
    d_(subDict.lookup<scalar>("d", sqrt(dimensions::specificHeatCapacity))),
    p_({a_*a_, b_, -2*a_*c_, -a_*d_, -c_*c_/3, -0.5*c_*d_, -0.2*d_*d_}),
    hsRef_(p_.integral(max(1 - constant::thermodynamic::Tstd/Tc_, small)))
{
    if (subDict.found("type"))
    {
        const word functionType = subDict.lookup<word>("type");

        if (functionType != Function1s::NSRDS14::typeName)
        {
            FatalIOErrorInFunction(subDict)
                << "Specific heat capacity function for NSRDS14 thermodynamics "
                << "should be of type " << Function1s::NSRDS14::typeName
                << ". Function type given is " << functionType << "."
                << exit(FatalIOError);
        }
    }
}


template<class EquationOfState>
Foam::NSRDS14Thermo<EquationOfState>::NSRDS14Thermo
(
    const word& name,
    const dictionary& dict
)
:
    NSRDS14Thermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::NSRDS14Thermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    writeEntry
    (
        os,
        "thermodynamics",
        dictionary::entries
        (
            "hf", hf_,
            "sf", sf_,
            "Tc", Tc_,
            "a", a_,
            "b", b_,
            "c", c_,
            "d", d_
        )
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const NSRDS14Thermo<EquationOfState>& ns
)
{
    ns.write(os);
    return os;
}


// ************************************************************************* //
