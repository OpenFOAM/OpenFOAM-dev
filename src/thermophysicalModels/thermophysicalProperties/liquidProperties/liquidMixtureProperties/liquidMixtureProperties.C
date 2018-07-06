/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
        if (dict.isDict(components_[i]))
        {
            properties_.set
            (
                i,
                liquidProperties::New(dict.subDict(components_[i]))
            );
        }
        else
        {
            properties_.set
            (
                i,
                liquidProperties::New(components_[i])
            );
        }
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

Foam::autoPtr<Foam::liquidMixtureProperties>
Foam::liquidMixtureProperties::New
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

Foam::scalar Foam::liquidMixtureProperties::Tc(const scalarField& X) const
{
    scalar vTc = 0.0;
    scalar vc = 0.0;

    forAll(properties_, i)
    {
        scalar x1 = X[i]*properties_[i].Vc();
        vc += x1;
        vTc += x1*properties_[i].Tc();
    }

    return vTc/vc;
}


Foam::scalar Foam::liquidMixtureProperties::Tpt(const scalarField& X) const
{
    scalar Tpt = 0.0;

    forAll(properties_, i)
    {
        Tpt += X[i]*properties_[i].Tt();
    }

    return Tpt;
}


Foam::scalar Foam::liquidMixtureProperties::pvInvert
(
    const scalar p,
    const scalarField& X
) const
{
    // Set upper and lower bounds
    scalar Thi = Tc(X);
    scalar Tlo = Tpt(X);

    // Check for critical and solid phase conditions
    if (p >= pv(p, Thi, X))
    {
        return Thi;
    }
    else if (p < pv(p, Tlo, X))
    {
        WarningInFunction
            << "Pressure below triple point pressure: "
            << "p = " << p << " < Pt = " << pv(p, Tlo, X) <<  nl << endl;
        return -1;
    }

    // Set initial guess
    scalar T = (Thi + Tlo)*0.5;

    while ((Thi - Tlo) > 1.0e-4)
    {
        if ((pv(p, T, X) - p) <= 0.0)
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


Foam::scalar Foam::liquidMixtureProperties::Tpc(const scalarField& X) const
{
    scalar Tpc = 0.0;

    forAll(properties_, i)
    {
        Tpc += X[i]*properties_[i].Tc();
    }

    return Tpc;
}


Foam::scalar Foam::liquidMixtureProperties::Ppc(const scalarField& X) const
{
    scalar Vc = 0.0;
    scalar Zc = 0.0;

    forAll(properties_, i)
    {
        Vc += X[i]*properties_[i].Vc();
        Zc += X[i]*properties_[i].Zc();
    }

    return RR*Zc*Tpc(X)/Vc;
}


Foam::scalar Foam::liquidMixtureProperties::omega(const scalarField& X) const
{
    scalar omega = 0.0;

    forAll(properties_, i)
    {
        omega += X[i]*properties_[i].omega();
    }

    return omega;
}


Foam::scalarField Foam::liquidMixtureProperties::Xs
(
    const scalar p,
    const scalar Tg,
    const scalar Tl,
    const scalarField& Xg,
    const scalarField& Xl
) const
{
    scalarField Xs(Xl.size());

    // Raoult's Law
    forAll(Xs, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), Tl);
        Xs[i] = properties_[i].pv(p, Ti)*Xl[i]/p;
    }

    return Xs;
}


Foam::scalar Foam::liquidMixtureProperties::W(const scalarField& X) const
{
    scalar W = 0.0;

    forAll(properties_, i)
    {
        W += X[i]*properties_[i].W();
    }

    return W;
}


Foam::scalarField Foam::liquidMixtureProperties::Y(const scalarField& X) const
{
    scalarField Y(X.size());
    scalar sumY = 0.0;

    forAll(Y, i)
    {
        Y[i] = X[i]*properties_[i].W();
        sumY += Y[i];
    }

    Y /= sumY;

    return Y;
}


Foam::scalarField Foam::liquidMixtureProperties::X(const scalarField& Y) const
{
    scalarField X(Y.size());
    scalar sumX = 0.0;

    forAll(X, i)
    {
        X[i] = Y[i]/properties_[i].W();
        sumX += X[i];
    }

    X /= sumX;

    return X;
}


Foam::scalar Foam::liquidMixtureProperties::rho
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    scalar sumY = 0.0;
    scalar v = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            scalar rho = properties_[i].rho(p, Ti);

            if (rho > small)
            {
                scalar Yi = X[i]*properties_[i].W();
                sumY += Yi;
                v += Yi/rho;
            }
        }
    }

    return sumY/v;
}


Foam::scalar Foam::liquidMixtureProperties::pv
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    scalar sumY = 0.0;
    scalar pv = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Yi = X[i]*properties_[i].W();
            sumY += Yi;

            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            pv += Yi*properties_[i].pv(p, Ti);
        }
    }

    return pv/sumY;
}


Foam::scalar Foam::liquidMixtureProperties::hl
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    scalar sumY = 0.0;
    scalar hl = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Yi = X[i]*properties_[i].W();
            sumY += Yi;

            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            hl += Yi*properties_[i].hl(p, Ti);
        }
    }

    return hl/sumY;
}


Foam::scalar Foam::liquidMixtureProperties::Cp
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    scalar sumY = 0.0;
    scalar Cp = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Yi = X[i]*properties_[i].W();
            sumY += Yi;

            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            Cp += Yi*properties_[i].Cp(p, Ti);
        }
    }

    return Cp/sumY;
}


Foam::scalar Foam::liquidMixtureProperties::sigma
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    // sigma is based on surface mole fractions
    // which are estimated from Raoult's Law
    scalar sigma = 0.0;
    scalarField Xs(X.size());
    scalar XsSum = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);
        scalar Pvs = properties_[i].pv(p, Ti);

        Xs[i] = X[i]*Pvs/p;
        XsSum += Xs[i];
    }

    Xs /= XsSum;

    forAll(properties_, i)
    {
        if (Xs[i] > small)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            sigma += Xs[i]*properties_[i].sigma(p, Ti);
        }
    }

    return sigma;
}


Foam::scalar Foam::liquidMixtureProperties::mu
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    scalar mu = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            mu += X[i]*log(properties_[i].mu(p, Ti));
        }
    }

    return exp(mu);
}


Foam::scalar Foam::liquidMixtureProperties::kappa
(
    const scalar p,
    const scalar T,
    const scalarField& X
) const
{
    // Calculate superficial volume fractions phii
    scalarField phii(X.size());
    scalar pSum = 0.0;

    forAll(properties_, i)
    {
        scalar Ti = min(TrMax*properties_[i].Tc(), T);

        scalar Vi = properties_[i].W()/properties_[i].rho(p, Ti);
        phii[i] = X[i]*Vi;
        pSum += phii[i];
    }

    phii /= pSum;

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
                    1.0/properties_[i].kappa(p, Ti)
                  + 1.0/properties_[j].kappa(p, Tj)
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
    const scalarField& X
) const
{
    // Blanc's law
    scalar Dinv = 0.0;

    forAll(properties_, i)
    {
        if (X[i] > small)
        {
            scalar Ti = min(TrMax*properties_[i].Tc(), T);
            Dinv += X[i]/properties_[i].D(p, Ti);
        }
    }

    return 1.0/Dinv;
}


// ************************************************************************* //
