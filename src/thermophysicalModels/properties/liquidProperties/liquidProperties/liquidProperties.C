/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "HashTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquidProperties, 0);
    defineRunTimeSelectionTable(liquidProperties,);
    defineRunTimeSelectionTable(liquidProperties, Istream);
    defineRunTimeSelectionTable(liquidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidProperties::liquidProperties
(
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
    W_(W),
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


Foam::liquidProperties::liquidProperties(Istream& is)
:
    W_(readScalar(is)),
    Tc_(readScalar(is)),
    Pc_(readScalar(is)),
    Vc_(readScalar(is)),
    Zc_(readScalar(is)),
    Tt_(readScalar(is)),
    Pt_(readScalar(is)),
    Tb_(readScalar(is)),
    dipm_(readScalar(is)),
    omega_(readScalar(is)),
    delta_(readScalar(is))
{}


Foam::liquidProperties::liquidProperties(const dictionary& dict)
:
    W_(readScalar(dict.lookup("W"))),
    Tc_(readScalar(dict.lookup("Tc"))),
    Pc_(readScalar(dict.lookup("Pc"))),
    Vc_(readScalar(dict.lookup("Vc"))),
    Zc_(readScalar(dict.lookup("Zc"))),
    Tt_(readScalar(dict.lookup("Tt"))),
    Pt_(readScalar(dict.lookup("Pt"))),
    Tb_(readScalar(dict.lookup("Tb"))),
    dipm_(readScalar(dict.lookup("dipm"))),
    omega_(readScalar(dict.lookup("omega"))),
    delta_(readScalar(dict.lookup("delta")))
{}


Foam::liquidProperties::liquidProperties(const liquidProperties& liq)
:
    W_(liq.W_),
    Tc_(liq.Tc_),
    Pc_(liq.Pc_),
    Vc_(liq.Vc_),
    Zc_(liq.Zc_),
    Tt_(liq.Tt_),
    Pt_(liq.Pt_),
    Tb_(liq.Tb_),
    dipm_(liq.dipm_),
    omega_(liq.omega_),
    delta_(liq.delta_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New(Istream& is)
{
    if (debug)
    {
        Info<< "liquidProperties::New(Istream&): "
            << "constructing liquidProperties" << endl;
    }

    const word liquidPropertiesType(is);
    const word coeffs(is);

    if (coeffs == "defaultCoeffs")
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(liquidPropertiesType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("liquidProperties::New(Istream&)")
                << "Unknown liquidProperties type "
                << liquidPropertiesType << nl << nl
                << "Valid liquidProperties types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquidProperties>(cstrIter()());
    }
    else if (coeffs == "coeffs")
    {
        return autoPtr<liquidProperties>(new liquidProperties(is));
    }
    else
    {
        FatalErrorIn("liquidProperties::New(Istream&)")
            << "liquidProperties type " << liquidPropertiesType
            << ", option " << coeffs << " given"
            << ", should be coeffs or defaultCoeffs"
            << abort(FatalError);

        return autoPtr<liquidProperties>(NULL);
    }
}


Foam::autoPtr<Foam::liquidProperties> Foam::liquidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "liquidProperties::New(const dictionary&):"
            << "constructing liquidProperties" << endl;
    }

    const word& liquidPropertiesTypeName = dict.dictName();

    const Switch defaultCoeffs(dict.lookup("defaultCoeffs"));

    if (defaultCoeffs)
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(liquidPropertiesTypeName);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "liquidProperties::New(const dictionary&)"
            )   << "Unknown liquidProperties type "
                << liquidPropertiesTypeName << nl << nl
                << "Valid liquidProperties types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquidProperties>(cstrIter()());
    }
    else
    {
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(liquidPropertiesTypeName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "liquidProperties::New(const dictionary&)"
            )   << "Unknown liquidProperties type "
                << liquidPropertiesTypeName << nl << nl
                << "Valid liquidProperties types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquidProperties>
        (
            cstrIter()(dict.subDict(liquidPropertiesTypeName + "Coeffs"))
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquidProperties::rho(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::rho(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::pv(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::pv(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::hl(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::hl(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::Cp(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::Cp(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::h(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::h(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::Cpg(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::Cpg(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::mu(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::mu(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::mug(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::mug(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::K(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::K(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::Kg(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::Kg(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::sigma(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::sigms(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::D(scalar p, scalar T) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::D(scalar, scalar) const"
    );
    return 0.0;
}


Foam::scalar Foam::liquidProperties::D(scalar p, scalar T, scalar Wb) const
{
    notImplemented
    (
        "Foam::scalar Foam::liquidProperties::D(scalar, scalar) const"
    );
    return 0.0;
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
            WarningIn
            (
                "Foam::scalar Foam::liquidProperties::pvInvert(scalar) const"
            )   << "Pressure below triple point pressure: "
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
        if ((pv(p, T) - p) <= 0.0)
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


void Foam::liquidProperties::writeData(Ostream& os) const
{

    os  << W_ << token::SPACE
        << Tc_ << token::SPACE
        << Pc_ << token::SPACE
        << Vc_ << token::SPACE
        << Zc_ << token::SPACE
        << Tt_ << token::SPACE
        << Pt_ << token::SPACE
        << Tb_ << token::SPACE
        << dipm_ << token::SPACE
        << omega_<< token::SPACE
        << delta_;
}


// ************************************************************************* //
