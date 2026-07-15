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

#include "janafThermo.H"
#include "units.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorInFunction
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ <= Tlow_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") <= Tlow(" << Tlow_ << ')'
            << exit(FatalError);
    }

    if (Tcommon_ > Thigh_)
    {
        FatalErrorInFunction
            << "Tcommon(" << Tcommon_ << ") > Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::janafThermo<EquationOfState>::janafThermo
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    EquationOfState(name, dict),
    Tlow_(subDict.lookup<scalar>("Tlow", dimensions::temperature)),
    Thigh_(subDict.lookup<scalar>("Thigh", dimensions::temperature)),
    Tcommon_(subDict.lookup<scalar>("Tcommon", dimensions::temperature)),
    highCpCoeffs_(subDict.lookup<coeffArray>("highCpCoeffs", units::none)),
    lowCpCoeffs_(subDict.lookup<coeffArray>("lowCpCoeffs", units::none))
{
    highCpCoeffs_ *= this->R();
    lowCpCoeffs_ *= this->R();

    checkInputData();
}


template<class EquationOfState>
Foam::janafThermo<EquationOfState>::janafThermo
(
    const word& name,
    const dictionary& dict
)
:
    janafThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janafThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    writeEntry
    (
        os,
        "thermodynamics",
        dictionary::entries
        (
            "Tlow", Tlow_,
            "Thigh", Thigh_,
            "Tcommon", Tcommon_,
            "highCpCoeffs", highCpCoeffs_/this->R(),
            "lowCpCoeffs", lowCpCoeffs_/this->R()
        )
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janafThermo<EquationOfState>& jt
)
{
    jt.write(os);
    return os;
}


// ************************************************************************* //
