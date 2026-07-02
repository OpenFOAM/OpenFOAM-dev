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

#include "NSRDSTransport.H"
#include "NSRDS0.H"
#include "NSRDS1.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::NSRDSTransport<Thermo>::NSRDSTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& muDict,
    const dictionary& kappaDict
)
:
    Thermo(name, dict),
    muA_(muDict.lookup<scalar>("a", units::none)),
    muB_(muDict.lookup<scalar>("b", units::none)),
    muC_(muDict.lookup<scalar>("c", units::none)),
    muD_(muDict.lookup<scalar>("d", units::none)),
    muE_(muDict.lookup<scalar>("e", units::none)),
    kappaCoeffs_("abcde", {dimTemperature, dimThermalConductivity}, kappaDict)
{
    if (muDict.found("type"))
    {
        const word functionType = muDict.lookup<word>("type");

        if (functionType != Function1s::NSRDS1::typeName)
        {
            FatalIOErrorInFunction(muDict)
                << "Viscosity function for NSRDS transport "
                << "should be of type " << Function1s::NSRDS1::typeName
                << ". Function type given is " << functionType << "."
                << exit(FatalIOError);
        }
    }

    if (kappaDict.found("type"))
    {
        const word functionType = kappaDict.lookup<word>("type");

        if (functionType != Function1s::NSRDS0::typeName)
        {
            FatalIOErrorInFunction(kappaDict)
                << "Thermal conductivity function for NSRDS transport "
                << "should be of type " << Function1s::NSRDS0::typeName
                << ". Function type given is " << functionType << "."
                << exit(FatalIOError);
        }
    }
}


template<class Thermo>
Foam::NSRDSTransport<Thermo>::NSRDSTransport
(
    const word& name,
    const dictionary& dict
)
:
    NSRDSTransport
    (
        name,
        dict,
        dict.subDict("transport").subDict("mu"),
        dict.subDict("transport").subDict("kappa")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::NSRDSTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    const delimitDictionary delimit(os, "transport");
    writeEntry
    (
        os,
        "mu",
        dictionary::entries
        (
            "a", muA_,
            "b", muB_,
            "c", muC_,
            "d", muD_,
            "e", muE_
        )
    );
    const delimitDictionary delimitKappa(os, "kappa");
    kappaCoeffs_.write("abcde", {dimTemperature, dimThermalConductivity}, os);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const NSRDSTransport<Thermo>& ns
)
{
    ns.write(os);
    return os;
}


// ************************************************************************* //
