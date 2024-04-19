/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "LangmuirHinshelwoodReactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const speciesTable& st,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    reactantNames_(dict.lookup("reactants")),
    a_(dict.lookupOrDefault<scalar>("a", 1)),
    A_(dict.lookup("A")),
    Ta_(dict.lookup("Ta")),
    beta_
    (
        dict.lookupOrDefault<FixedList<scalar, 3>>
        (
            "beta",
            FixedList<scalar, 3>({0, 0, 0})
        )
    ),
    m_
    (
        dict.lookupOrDefault<FixedList<scalar, 3>>
        (
            "m",
            FixedList<scalar, 3>({2, 1, 1})
        )
    )
{
    forAll(reactantNames_, i)
    {
        r_[i] = st[reactantNames_[i]];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::LangmuirHinshelwoodReactionRate::preEvaluate() const
{}


inline void Foam::LangmuirHinshelwoodReactionRate::postEvaluate() const
{}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar c1m = pow(c[r_[0]], m_[1]);
    const scalar c2m = pow(c[r_[1]], m_[2]);

    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-Ta_[0]/T);
    const scalar k1 = A_[1]*pow(T, beta_[1])*exp(-Ta_[1]/T);
    const scalar k2 = A_[2]*pow(T, beta_[2])*exp(-Ta_[2]/T);

    const scalar b = a_ + k1*c1m + k2*c2m;

    return k0/pow(b, m_[0]);
}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar c1m = pow(c[r_[0]], m_[1]);
    const scalar c2m = pow(c[r_[1]], m_[2]);

    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-Ta_[0]/T);
    const scalar k1 = A_[1]*pow(T, beta_[1])*exp(-Ta_[1]/T);
    const scalar k2 = A_[2]*pow(T, beta_[2])*exp(-Ta_[2]/T);

    const scalar dk0dT = k0/T*(beta_[0] + Ta_[0]/T);
    const scalar dk1dT = k1/T*(beta_[1] + Ta_[1]/T);
    const scalar dk2dT = k2/T*(beta_[2] + Ta_[2]/T);

    const scalar b = a_ + k1*c1m + k2*c2m;
    const scalar dbdT = dk1dT*c1m + dk2dT*c2m;

    return (dk0dT - k0*m_[0]*dbdT/b)/pow(b, m_[0]);
}


inline bool Foam::LangmuirHinshelwoodReactionRate::hasDdc() const
{
    return true;
}


inline void Foam::LangmuirHinshelwoodReactionRate::ddc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& ddc
) const
{
    const scalar c1m = pow(c[r_[0]], m_[1]);
    const scalar c2m = pow(c[r_[1]], m_[2]);

    const scalar dc1mdc1 = m_[1]*c1m/c[r_[0]];
    const scalar dc2mdc2 = m_[2]*c2m/c[r_[1]];

    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-Ta_[0]/T);
    const scalar k1 = A_[1]*pow(T, beta_[1])*exp(-Ta_[1]/T);
    const scalar k2 = A_[2]*pow(T, beta_[2])*exp(-Ta_[2]/T);

    const scalar b = a_ + k1*c1m + k2*c2m;
    const scalar dbdc1 = k1*dc1mdc1;
    const scalar dbdc2 = k2*dc2mdc2;

    const scalar k = k0/pow(b, m_[0]);
    const scalar dkdb = - k*m_[0]/b;

    ddc = 0;
    ddc[r_[0]] = dkdb*dbdc1;
    ddc[r_[1]] = dkdb*dbdc2;
}


inline void Foam::LangmuirHinshelwoodReactionRate::write(Ostream& os) const
{
    writeEntry(os, "reactants", reactantNames_);
    writeEntry(os, "a", a_);
    writeEntry(os, "A", A_);
    writeEntry(os, "Ta", Ta_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "m", m_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LangmuirHinshelwoodReactionRate& lhrr
)
{
    lhrr.write(os);
    return os;
}


// ************************************************************************* //
