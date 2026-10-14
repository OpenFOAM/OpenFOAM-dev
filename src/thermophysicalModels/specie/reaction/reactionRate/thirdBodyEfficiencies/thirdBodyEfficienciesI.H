/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "thirdBodyEfficiencies.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const List<Tuple2<word, scalar>>& efficiencies
)
:
    species_(species),
    efficiencies_(species.size(), 1.0)
{
    nonUnitySpecies_.setSize(efficiencies.size());
    nonUnityEfficienciesM1_.setSize(efficiencies.size());
    label n = 0;

    forAll(efficiencies, i)
    {
        if (efficiencies[i].second() != 1)
        {
            efficiencies_[species[efficiencies[i].first()]] =
                efficiencies[i].second();

            nonUnitySpecies_[n] = species[efficiencies[i].first()];
            nonUnityEfficienciesM1_[n] = efficiencies[i].second() - 1;
            n++;
        }
    }

    nonUnitySpecies_.setSize(n);
    nonUnityEfficienciesM1_.setSize(n);
}


inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const dictionary& dict
)
:
    thirdBodyEfficiencies
    (
        species,
        dict.lookup<List<Tuple2<word, scalar>>>("coeffs")
    )
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M
(
    const reaction::state& s
) const
{
    scalar M = s.C;
    forAll(nonUnitySpecies_, i)
    {
        M += nonUnityEfficienciesM1_[i]*s.c[nonUnitySpecies_[i]];
    }

    return M;
}


inline const Foam::scalarField&
Foam::thirdBodyEfficiencies::dMdc() const
{
    return efficiencies_;
}


inline void Foam::thirdBodyEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> efficiencies(nonUnitySpecies_.size());
    forAll(nonUnitySpecies_, i)
    {
        efficiencies[i].first() = species_[nonUnitySpecies_[i]];
        efficiencies[i].second() = nonUnityEfficienciesM1_[i] + 1;
    }

    writeEntry(os, "coeffs", efficiencies);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thirdBodyEfficiencies& tbes
)
{
    tbes.write(os);
    return os;
}


// ************************************************************************* //
