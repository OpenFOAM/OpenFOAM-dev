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

#include "FixedList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const scalar A[],
    const scalar Ta[],
    const label co,
    const label c3h6,
    const label no
)
:
    co_(co),
    c3h6_(c3h6),
    no_(no)
{
    for (int i=0; i<n_; i++)
    {
        A_[i] = A[i];
        Ta_[i] = Ta[i];
    }
}


inline Foam::LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const speciesTable& st,
    const dictionary& dict
)
:
    co_(st["CO"]),
    c3h6_(st["C3H6"]),
    no_(st["NO"])
{
    // read (A, Ta) pairs
    FixedList<Tuple2<scalar, scalar>, n_> coeffs(dict.lookup("coeffs"));

    forAll(coeffs, i)
    {
        A_[i] = coeffs[i].first();
        Ta_[i] = coeffs[i].second();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return A_[0]*exp(-Ta_[0]/T)/
    (
        T
       *sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
    );
}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    scalar den =
    (
        T
       *sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
    );
    scalar rate = A_[0]*exp(-Ta_[0]/T)/ den;

    scalar derivDen =
    (
        sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
      + 2*T*(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *
        (
            A_[1]*exp(-Ta_[1]/T)*c[co_]*Ta_[1]/sqr(T)
          + A_[2]*exp(-Ta_[2]/T)*c[c3h6_]*Ta_[2]/sqr(T)
        )
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
      + T
       *sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(A_[3]*exp(-Ta_[3]/T)*Ta_[3]*sqr(c[co_])*sqr(c[c3h6_])/sqr(T))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
      + T
       *sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(A_[4]*exp(-Ta_[4]/T)*Ta_[4]*pow(c[no_], 0.7))/sqr(T)
    );

    return rate*(Ta_[0]/sqr(T) - derivDen/den);
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
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0;
}


inline void Foam::LangmuirHinshelwoodReactionRate::write(Ostream& os) const
{
    FixedList<Tuple2<scalar, scalar>, n_> coeffs;

    forAll(coeffs, i)
    {
        coeffs[i].first() = A_[i];
        coeffs[i].second() = Ta_[i];
    }

    os.writeKeyword("coeffs") << coeffs << nl;
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
