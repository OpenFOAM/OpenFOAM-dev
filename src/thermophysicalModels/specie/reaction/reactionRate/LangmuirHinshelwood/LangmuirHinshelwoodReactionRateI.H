/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "FixedList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const speciesTable& st,
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

    const scalar TaByT0 = Ta_[0]/T;
    const scalar TaByT1 = Ta_[1]/T;
    const scalar TaByT2 = Ta_[2]/T;

    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-TaByT0);
    const scalar k1 = A_[1]*pow(T, beta_[1])*exp(-TaByT1);
    const scalar k2 = A_[2]*pow(T, beta_[2])*exp(-TaByT2);

    return k0/pow(a_ + k1*c1m + k2*c2m, m_[0]);
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

    const scalar TaByT0 = Ta_[0]/T;
    const scalar TaByT1 = Ta_[1]/T;
    const scalar TaByT2 = Ta_[2]/T;

    const scalar k0 = A_[0]*pow(T, beta_[0])*exp(-TaByT0);
    const scalar k1 = A_[1]*pow(T, beta_[1])*exp(-TaByT1);
    const scalar k2 = A_[2]*pow(T, beta_[2])*exp(-TaByT2);

    return
       (
           (beta_[0] + TaByT0)*k0
         - m_[0]*k0*((beta_[1] + TaByT1)*k1*c1m + (beta_[2] + TaByT2)*k2*c2m)
          /(a_ + k1*c1m + k2*c2m)
       )/(pow(a_ + k1*c1m + k2*c2m, m_[0])*T);
}


inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::LangmuirHinshelwoodReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}


inline void Foam::LangmuirHinshelwoodReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
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
