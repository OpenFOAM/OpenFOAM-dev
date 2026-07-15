/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "dictionary.H"
#include <cfenv>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict

)
:
    Specie(name, dict),
    Tc_(subDict.lookup<scalar>("Tc", dimensions::temperature)),
    Vc_(subDict.lookup<scalar>("Vc", dimensions::volume/dimensions::moles)),
    Zc_(1.0),
    Pc_(subDict.lookup<scalar>("Pc", dimensions::pressure)),
    omega_(subDict.lookup<scalar>("omega", dimless))
{
    Zc_ = Pc_*Vc_/(constant::thermodynamic::RR*Tc_);
}


template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const word& name,
    const dictionary& dict
)
:
    PengRobinsonGas(name, dict, dict.subDict("equationOfState"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinsonGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    writeEntry
    (
        os,
        "equationOfState",
        dictionary::entries("Tc", Tc_, "Vc", Vc_, "Pc", Pc_, "omega", omega_)
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonGas<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
