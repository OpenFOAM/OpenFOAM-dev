/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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

#include "volFields.H"
#include "mathematicalConstants.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::mathematical;
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar Foam::fluxLimitedLangmuirHinshelwoodReactionRate::Av
(
    const label li
) const
{
    if (AvUniform_)
    {
        return Av_;
    }
    else
    {
        return AvField_[li];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::fluxLimitedLangmuirHinshelwoodReactionRate::
fluxLimitedLangmuirHinshelwoodReactionRate
(
    const speciesTable& st,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    additionalAdsorbableSpecieNames_
    (
        dict.lookupOrDefault
        (
            "additionalAdsorbableSpecies",
            wordList()
        )
    ),
    a_(dict.lookupOrDefault("a", scalar(1))),
    A_(dict.lookup("A")),
    Ta_(dict.lookup("Ta")),
    limited_(dict.found("s")),
    AvUniform_
    (
        limited_
      ? dict.lookup("Av")[0].isNumber()
      : true
    ),
    Av_
    (
        limited_ && AvUniform_
      ? dict.lookup<scalar>("Av")
      : 0
    ),
    AvName_
    (
        limited_ && !AvUniform_
      ? dict.lookup("Av")
      : word::null
    ),
    AvField_
    (
        limited_ && !AvUniform_
      ? ob.lookupObject<volScalarField::Internal>(AvName_)
      : volScalarField::null()
    )
{
    List<specieCoeffs> lhs;
    List<specieCoeffs> rhs;
    specieCoeffs::setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        st,
        lhs,
        rhs
    );

    nReactants_ = lhs.size();

    ra_.setSize(nReactants_ + additionalAdsorbableSpecieNames_.size());
    nu_.setSize(nReactants_);
    exp_.setSize(nReactants_);

    forAll(lhs, i)
    {
        ra_[i] = lhs[i].index;
        nu_[i] = lhs[i].stoichCoeff;
        exp_[i] = lhs[i].exponent;
    }

    forAll(additionalAdsorbableSpecieNames_, i)
    {
        ra_[nReactants_ + i] = st[additionalAdsorbableSpecieNames_[i]];
    }

    beta_ =
        dict.lookupOrDefault
        (
            "beta",
            scalarList
            (
                1 + nReactants_ + additionalAdsorbableSpecieNames_.size(),
                scalar(0)
            )
        );

    m_ =
        dict.lookupOrDefault
        (
            "m",
            scalarList
            (
                1 + nReactants_ + additionalAdsorbableSpecieNames_.size(),
                scalar(1)
            )
        );

    if (!dict.found("m"))
    {
        m_[0] = 2;
    }

    s_ = dict.lookupOrDefault("s", scalarList(nReactants_, scalar(0)));

    W_ =
        limited_
      ? scalarList(dict.lookup("W"))
      : scalarList(nReactants_, scalar(0));

    const scalar nCoeffs =
        1 + nReactants_ + additionalAdsorbableSpecieNames_.size();

    if (A_.size() != nCoeffs)
    {
        FatalIOErrorInFunction(dict)
           << "Number of A coefficients != " << nCoeffs
           << exit(FatalIOError);
    }

    if (Ta_.size() != nCoeffs)
    {
        FatalIOErrorInFunction(dict)
           << "Number of Ta coefficients != " << nCoeffs
           << exit(FatalIOError);
    }

    if (beta_.size() != nCoeffs)
    {
        FatalIOErrorInFunction(dict)
           << "Number of beta coefficients != " << nCoeffs
           << exit(FatalIOError);
    }

    if (m_.size() != nCoeffs)
    {
        FatalIOErrorInFunction(dict)
           << "Number of m coefficients != " << nCoeffs
           << exit(FatalIOError);
    }

    if (s_.size() != nReactants_)
    {
        FatalIOErrorInFunction(dict)
           << "Number of s coefficients != "
           << nReactants_ << exit(FatalIOError);
    }

    if (W_.size() != nReactants_)
    {
        FatalIOErrorInFunction(dict)
           << "Number of W coefficients != "
           << nReactants_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void
Foam::fluxLimitedLangmuirHinshelwoodReactionRate::preEvaluate() const
{}


inline void
Foam::fluxLimitedLangmuirHinshelwoodReactionRate::postEvaluate() const
{}


inline Foam::scalar
Foam::fluxLimitedLangmuirHinshelwoodReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    scalar sumkc = 0;
    forAll(ra_, i)
    {
        const label ip1 = i + 1;

        const scalar kc =
            A_[ip1]*pow(T, beta_[ip1])*exp(-Ta_[ip1]/T)
           *pow(c[ra_[i]], m_[ip1]);

        sumkc += kc;
    }

    const scalar TaByT0 = Ta_[0]/T;
    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-TaByT0);

    scalar r = k0/pow(a_ + sumkc, m_[0]);

    if (limited_)
    {
        scalar rc = 1;
        forAll(ra_, i)
        {
            rc *= pow(c[ra_[i]], exp_[i]);
        }

        const scalar Av = this->Av(li);

        if (rc > rootSmall)
        {
            forAll(s_, i)
            {
                r = min
                (
                    r,
                    (Av*s_[i]*c[ra_[i]]/(nu_[i]*rc))*sqrt(RR*T/(twoPi*W_[i]))
                );
            }
        }
    }

    return r;
}


inline Foam::scalar Foam::fluxLimitedLangmuirHinshelwoodReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    scalar sumkc = 0;
    scalar sumBetaKc = 0;
    forAll(ra_, i)
    {
        const label ip1 = i + 1;

        const scalar kc =
            A_[ip1]*pow(T, beta_[ip1])*exp(-Ta_[ip1]/T)
           *pow(c[ra_[i]], m_[ip1]);

        sumkc += kc;
        sumBetaKc += (beta_[ip1] + Ta_[ip1]/T)*kc;
    }

    const scalar TaByT0 = Ta_[0]/T;
    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-TaByT0);

    scalar ddT =
       ((beta_[0] + TaByT0)*k0 - m_[0]*k0*sumBetaKc/(a_ + sumkc))
      /(pow(a_ + sumkc, m_[0])*T);

    if (limited_)
    {
        scalar rc = 1;
        forAll(ra_, i)
        {
            rc *= pow(c[ra_[i]], exp_[i]);
        }

        scalar r = k0/pow(a_ + sumkc, m_[0]);

        if (rc > rootSmall)
        {
            label l = -1;

            const scalar Av = this->Av(li);

            forAll(s_, i)
            {
                const scalar rl =
                    (Av*s_[i]*c[ra_[i]]/(nu_[i]*rc))*sqrt(RR*T/(twoPi*W_[i]));

                if (rl < r)
                {
                    l = i;
                    r = rl;
                }
            }

            if (l != -1)
            {
                ddT =
                    (Av*s_[l]*c[ra_[l]]/(nu_[l]*rc))
                   *0.5*sqrt(RR/(twoPi*W_[l]*T));
            }
        }
    }

    return ddT;
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::fluxLimitedLangmuirHinshelwoodReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::fluxLimitedLangmuirHinshelwoodReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::fluxLimitedLangmuirHinshelwoodReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline void Foam::fluxLimitedLangmuirHinshelwoodReactionRate::write
(
    Ostream& os
) const
{
    writeEntry
    (
        os,
        "additionalAdsorbableSpecies",
        additionalAdsorbableSpecieNames_
    );
    writeEntry(os, "a", a_);
    writeEntry(os, "A", A_);
    writeEntry(os, "Ta", Ta_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "m", m_);
    if (limited_)
    {
        writeEntry(os, "s", s_);
        writeEntry(os, "W", W_);
    }
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fluxLimitedLangmuirHinshelwoodReactionRate& lhrr
)
{
    lhrr.write(os);
    return os;
}


// ************************************************************************* //
