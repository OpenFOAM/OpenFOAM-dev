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

#include "ePolynomialThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
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
    CvCoeffs_
    (
        subDict.lookup<FixedLaurentPolynomial<scalar, 0, PolySize>>
        (
            "CvCoeffs<" + Foam::name(PolySize) + '>',
            Function1s::unitSets({dimTemperature, dimSpecificHeatCapacity})
        )
    ),
    esRef_(CvCoeffs_.integral(constant::thermodynamic::Tstd)),
    sRef_(CvCoeffs_.byX().integral(constant::thermodynamic::Tstd))
{}


template<class EquationOfState, int PolySize>
Foam::ePolynomialThermo<EquationOfState, PolySize>::ePolynomialThermo
(
    const word& name,
    const dictionary& dict
)
:
    ePolynomialThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
void Foam::ePolynomialThermo<EquationOfState, PolySize>::write
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
            word("CvCoeffs<" + Foam::name(PolySize) + '>'), CvCoeffs_
        )
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ePolynomialThermo<EquationOfState, PolySize>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
