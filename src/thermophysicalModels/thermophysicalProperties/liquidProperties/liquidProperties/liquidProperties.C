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

#include "liquidProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidProperties, 0);
    defineRunTimeSelectionTable(liquidProperties,);
    defineRunTimeSelectionTable(liquidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidProperties::liquidProperties
(
    const word& name,
    scalar W,
    scalar Tc,
    scalar Pc,
    scalar Vc,
    scalar Zc,
    scalar Tt,
    scalar Pt,
    scalar Tb,
    scalar dipm,
    scalar omega,
    scalar delta
)
:
    thermophysicalProperties(W),
    name_(name),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    Zc_(Zc),
    Tt_(Tt),
    Pt_(Pt),
    Tb_(Tb),
    dipm_(dipm),
    omega_(omega),
    delta_(delta)
{}


Foam::liquidProperties::liquidProperties(const dictionary& dict)
:
    thermophysicalProperties(dict),
    name_(dict.dictName()),
    Tc_(dict.lookup<scalar>("Tc")),
    Pc_(dict.lookup<scalar>("Pc")),
    Vc_(dict.lookup<scalar>("Vc")),
    Zc_(dict.lookup<scalar>("Zc")),
    Tt_(dict.lookup<scalar>("Tt")),
    Pt_(dict.lookup<scalar>("Pt")),
    Tb_(dict.lookup<scalar>("Tb")),
    dipm_(dict.lookup<scalar>("dipm")),
    omega_(dict.lookup<scalar>("omega")),
    delta_(dict.lookup<scalar>("delta"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::liquidProperties::name() const
{
    return name_;
}


Foam::scalar Foam::liquidProperties::S(scalar p, scalar T) const
{
    NotImplemented;
    return 0;
}


Foam::scalar Foam::liquidProperties::pvInvert(scalar p) const
{
    // Check for critical and solid phase conditions
    if (p >= Pc_)
    {
        return Tc_;
    }
    else if (p < Pt_)
    {
        if (debug)
        {
            WarningInFunction
                << "Pressure below triple point pressure: "
                << "p = " << p << " < Pt = " << Pt_ <<  nl << endl;
        }
        return -1;
    }

    // Set initial upper and lower bounds
    scalar Thi = Tc_;
    scalar Tlo = Tt_;

    // Initialise T as boiling temperature under normal conditions
    scalar T = Tb_;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T) - p) <= 0)
        {
            Tlo = T;
        }
        else
        {
            Thi = T;
        }

        T = (Thi + Tlo)*0.5;
    }

    return T;
}


void Foam::liquidProperties::readIfPresent(const dictionary &dict)
{
    thermophysicalProperties::readIfPresent(dict);
    dict.readIfPresent("Tc", Tc_);
    dict.readIfPresent("Pc", Pc_);
    dict.readIfPresent("Vc", Vc_);
    dict.readIfPresent("Zc", Zc_);
    dict.readIfPresent("Tt", Tt_);
    dict.readIfPresent("Pt", Pt_);
    dict.readIfPresent("Tb", Tb_);
    dict.readIfPresent("dipm", dipm_);
    dict.readIfPresent("omega", omega_);
    dict.readIfPresent("delta", delta_);
}


void Foam::liquidProperties::write(Ostream& os) const
{
    thermophysicalProperties::write(os);
    writeEntry(os, "Tc", Tc_);
    writeEntry(os, "Pc", Pc_);
    writeEntry(os, "Vc", Vc_);
    writeEntry(os, "Zc", Zc_);
    writeEntry(os, "Tt", Tt_);
    writeEntry(os, "Pt", Pt_);
    writeEntry(os, "Tb", Tb_);
    writeEntry(os, "dipm", dipm_);
    writeEntry(os, "omega", omega_);
    writeEntry(os, "delta", delta_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const liquidProperties& l)
{
    l.write(os);
    return os;
}


// ************************************************************************* //
