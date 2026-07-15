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

#include "NSRDSThermo.H"
#include "NSRDS0.H"
#include "thermodynamicConstants.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::NSRDSThermo<EquationOfState>::NSRDSThermo
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
    CpCoeffs_
    (
        "abcde",
        {dimensions::temperature, dimensions::specificHeatCapacity},
        subDict
    ),
    hsRef_(CpCoeffs_.integral(constant::thermodynamic::Tstd)),
    sRef_(CpCoeffs_.byX().integral(constant::thermodynamic::Tstd))
{
    if (subDict.found("type"))
    {
        const word functionType = subDict.lookup<word>("type");

        if (functionType != Function1s::NSRDS0::typeName)
        {
            FatalIOErrorInFunction(subDict)
                << "Specific heat capacity function for NSRDS thermodynamics "
                << "should be of type " << Function1s::NSRDS0::typeName
                << ". Function type given is " << functionType << "."
                << exit(FatalIOError);
        }
    }
}


template<class EquationOfState>
Foam::NSRDSThermo<EquationOfState>::NSRDSThermo
(
    const word& name,
    const dictionary& dict
)
:
    NSRDSThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::NSRDSThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    const delimitDictionary delimit(os, "thermodynamics");
    writeEntry(os, "hf", hf_);
    writeEntry(os, "sf", sf_);
    CpCoeffs_.write
    (
        "abcde",
        {dimensions::temperature, dimensions::specificHeatCapacity},
        os
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const NSRDSThermo<EquationOfState>& ns
)
{
    ns.write(os);
    return os;
}


// ************************************************************************* //
