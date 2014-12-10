/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "liquidMixtureProperties.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::liquidMixtureProperties::TrMax = 0.999;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquidMixtureProperties::liquidMixtureProperties
(
    const dictionary& dict
)
:
    components_(),
    properties_()
{
    components_ = dict.toc();
    properties_.setSize(components_.size());

    forAll(components_, i)
    {
        properties_.set
        (
            i,
            liquidProperties::New(dict.subDict(components_[i]))
        );
    }
}


Foam::liquidMixtureProperties::liquidMixtureProperties
(
    const liquidMixtureProperties& lm
)
:
    components_(lm.components_),
    properties_(lm.properties_.size())
{
    forAll(properties_, i)
    {
        properties_.set(i, lm.properties_(i)->clone());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquidMixtureProperties> Foam::liquidMixtureProperties::New
(
    const dictionary& thermophysicalProperties
)
{
    return autoPtr<liquidMixtureProperties>
    (
        new liquidMixtureProperties(thermophysicalProperties)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquidMixtureProperties::Tc(const scalarField& x) const
{
    scalar vTc = 0.0;
    scalar vc = 0.0;

    forAll(properties_, i)
    {
        scalar x1 = x[i]*properties_[i].Vc();
        vc += x1;
        vTc += x1*properties_[i].Tc();
    }

    return vTc/vc;
}


Foam::scalar Foam::liquidMixtureProperties::Tpt(const scalarField& x) const
{
    scalar Tpt = 0.0;
    forAll(properties_, i)
    {
        Tpt += x[i]*properties_[i].Tt();
    }

    return Tpt;
}


Foam::scalar Foam::liquidMixtureProperties::pvInvert
(
    const scalar p,
    const scalarField& x
) const
{
    // Set upper and lower bounds
    scalar Thi = Tc(x);
    scalar Tlo = Tpt(x);

    // Check for critical and solid phase conditions
    if (p >= pv(p, Thi, x))
    {
        return Thi;
    }
    else if (p < pv(p, Tlo, x))
    {
        WarningIn
        (
            "Foam::scalar Foam::liquidMixtureProperties::pvInvert"
            "("
            "    const scalar,"
            "    const scalarField&"
            ") const"
        )   << "Pressure below triple point pressure: "
            << "p = " << p << " < Pt = " << pv(p, Tlo, x) <<  nl << endl;
        return -1;
    }

    // Set initial guess
    scalar T = (Thi + Tlo)*0.5;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T, x) - p) <= 0.0)
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


Foam::scalar Foam::liquidMixtureProperties::Tpc(const scalarField& x) const
{
    scalar Tpc = 0.0;
    forAll(properties_, i)
    {
        Tpc += x[i]*properties_[i].Tc();
    }

    return Tpc;
}


Foam::scalar Foam::liquidMixtureProperties::Ppc(const scalarField& x) const
{
    scalar Vc = 0.0;
    scalar Zc = 0.0;
    forAll(properties_, i)
    {
        Vc += x[i]*properties_[i].Vc();
        Zc += x[i]*properties_[i].Zc();
    }

    return specie::RR*Zc*Tpc(x)/Vc;
}


Foam::scalar Foam::liquidMixtureProperties::omega(const scalarField& x) const
{
    scalar omega = 0.0;
    forAll(properties_, i)
    {
        omega += x[i]*properties_[i].omega();
    }

    return omega;
}


Foam::scalarField Foam::liquidMixtureProperties::Xs
(
    const scalar p,
    const scalar Tg,
    const scalar Tl,
    const scalarField& xg,
    const scalarField& xl
) const
{
    scalarField xs(xl.size(), 0.0);

    // Raoult's Law
    forAll(xs, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), Tl);
        xs[i] = properties_[i].pv(p, Ti)*xl[i]/p;
    }
    return xs;
}


Foam::scalar Foam::liquidMixtureProperties::W(const scalarField& x) const
{
    scalar W = 0.0;
    forAll(properties_, i)
    {
        W += x[i]*properties_[i].W();
    }

    return W;
}


Foam::scalarField Foam::liquidMixtureProperties::Y(const scalarField& X) const
{
    scalarField Y(X/W(X));

    forAll(Y, i)
    {
        Y[i] *= properties_[i].W();
    }

    return Y;
}


Foam::scalarField Foam::liquidMixtureProperties::X(const scalarField& Y) const
{
    scalarField X(Y.size());
    scalar Winv = 0.0;
    forAll(X, i)
    {
        Winv += Y[i]/properties_[i].W();
        X[i] = Y[i]/properties_[i].W();
    }

    tmp<scalarField> tfld = X/Winv;
    return tfld();
}


Foam::scalar Foam::liquidMixtureProperties::rho
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar v = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            scalar rho = SMALL + properties_[i].rho(p, Ti);
            v += x[i]*properties_[i].W()/rho;
        }
    }

    return W(x)/v;
}


Foam::scalar Foam::liquidMixtureProperties::pv
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar pv = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            pv += x[i]*properties_[i].pv(p, Ti)*properties_[i].W();
        }
    }

    return pv/W(x);
}


Foam::scalar Foam::liquidMixtureProperties::hl
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar hl = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            hl += x[i]*properties_[i].hl(p, Ti)*properties_[i].W();
        }
    }

    return hl/W(x);
}


Foam::scalar Foam::liquidMixtureProperties::Cp
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar Cp = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            Cp += x[i]*properties_[i].Cp(p, Ti)*properties_[i].W();
        }
    }

    return Cp/W(x);
}


Foam::scalar Foam::liquidMixtureProperties::sigma
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // sigma is based on surface mole fractions
    // which is estimated from Raoult's Law
    scalar sigma = 0.0;
    scalarField Xs(x.size(), 0.0);
    scalar XsSum = 0.0;
    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);
        scalar Pvs = properties_[i].pv(p, Ti);
        scalar xs = x[i]*Pvs/p;
        XsSum += xs;
        Xs[i] = xs;
    }

    forAll(properties_, i)
    {
        if (Xs[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            sigma += (Xs[i]/XsSum)*properties_[i].sigma(p, Ti);
        }
    }

    return sigma;
}


Foam::scalar Foam::liquidMixtureProperties::mu
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    scalar mu = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            mu += x[i]*log(properties_[i].mu(p, Ti));
        }
    }

    return exp(mu);
}


Foam::scalar Foam::liquidMixtureProperties::K
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // calculate superficial volume fractions phii
    scalarField phii(x.size(), 0.0);
    scalar pSum = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        scalar Vi = properties_[i].W()/properties_[i].rho(p, Ti);
        phii[i] = x[i]*Vi;
        pSum += phii[i];
    }

    forAll(phii, i)
    {
        phii[i] /= pSum;
    }

    scalar K = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        forAll(properties_, j)
        {
            scalar Tj = min(TrMax*properties_[j].Tc(), T);

            scalar Kij =
                2.0
               /(
                    1.0/properties_[i].K(p, Ti)
                  + 1.0/properties_[j].K(p, Tj)
                );
            K += phii[i]*phii[j]*Kij;
        }
    }

    return K;
}


Foam::scalar Foam::liquidMixtureProperties::D
(
    const scalar p,
    const scalar T,
    const scalarField& x
) const
{
    // Blanc's law
    scalar Dinv = 0.0;

    forAll(properties_, i)
    {
        if (x[i] > SMALL)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            Dinv += x[i]/properties_[i].D(p, Ti);
        }
    }

    return 1.0/Dinv;
}


// ************************************************************************* //
