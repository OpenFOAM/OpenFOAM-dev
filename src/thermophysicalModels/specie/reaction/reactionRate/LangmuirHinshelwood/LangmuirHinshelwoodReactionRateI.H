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
    const speciesTable& st,
    const dictionary& dict
)
:
    species_(st),
    reactantNames_(dict.lookup("reactants")),
    a_(dict.lookupOrDefault<scalar>("a", 1)),
    A_(dict.lookup("A")),
    Ta_(dict.lookup("Ta")),
    m_
    (
        dict.lookupOrDefault<FixedList<scalar, 3>>
        (
            "m",
            FixedList<scalar, 3>({1, 1, 2})
        )
    )
{
    forAll(reactantNames_, i)
    {
        r_[i] = st[reactantNames_[i]];
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
    const scalar c0m = pow(c[r_[0]], m_[0]);
    const scalar c1m = pow(c[r_[1]], m_[1]);

    return A_[0]*exp(-Ta_[0]/T)/
    (
        T*pow(a_ + A_[1]*exp(-Ta_[1]/T)*c0m + A_[2]*exp(-Ta_[2]/T)*c1m, m_[2])
    );
}


inline Foam::scalar Foam::LangmuirHinshelwoodReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    const scalar c0m = pow(c[r_[0]], m_[0]);
    const scalar c1m = pow(c[r_[1]], m_[1]);

    const scalar dr0 = A_[1]*exp(-Ta_[1]/T)*c0m;
    const scalar dr1 = A_[2]*exp(-Ta_[2]/T)*c1m;

    const scalar denByT = pow(a_ + dr0 + dr1, m_[2]);
    const scalar den = T*denByT;
    const scalar rate = A_[0]*exp(-Ta_[0]/T)/den;

    const scalar derivDen =
    (
        denByT
      + m_[2]*T*(a_ + dr0 + dr1)*(dr0*Ta_[1]/sqr(T) + dr1*Ta_[2]/sqr(T))
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
    os.writeKeyword("reactants")
        << reactantNames_ << token::END_STATEMENT << nl;
    os.writeKeyword("a") << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
    os.writeKeyword("m") << m_ << token::END_STATEMENT << nl;
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
