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
    const scalar defaultEfficiency,
    const List<Tuple2<word, scalar>>& efficiencies
)
:
    species_(species),
    defaultEfficiency_(defaultEfficiency),
    efficiencies_(species.size(), defaultEfficiency_)
{
    thirdBodySpecies_.setSize(efficiencies.size());
    efficienciesMdefault_.setSize(efficiencies.size());
    label n = 0;

    forAll(efficiencies, i)
    {
        if (efficiencies[i].second() != defaultEfficiency_)
        {
            efficiencies_[species[efficiencies[i].first()]] =
                efficiencies[i].second();

            thirdBodySpecies_[n] = species[efficiencies[i].first()];
            efficienciesMdefault_[n] =
                efficiencies[i].second() - defaultEfficiency_;
            n++;
        }
    }

    thirdBodySpecies_.setSize(n);
    efficienciesMdefault_.setSize(n);
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
        dict.lookup<scalar>("defaultEfficiency"),
        dict.lookupOrDefault<List<Tuple2<word, scalar>>>
        (
            "efficiencies",
            List<Tuple2<word, scalar>>()
        )
    )
{}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M
(
    const reaction::state& s
) const
{
    scalar M = defaultEfficiency_*s.C;
    forAll(thirdBodySpecies_, i)
    {
        M += efficienciesMdefault_[i]*s.c[thirdBodySpecies_[i]];
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
    writeEntry(os, "defaultEfficiency", defaultEfficiency_);

    if (thirdBodySpecies_.size())
    {
        List<Tuple2<word, scalar>> efficiencies(thirdBodySpecies_.size());

        forAll(thirdBodySpecies_, i)
        {
            efficiencies[i].first() = species_[thirdBodySpecies_[i]];
            efficiencies[i].second() =
                efficienciesMdefault_[i] + defaultEfficiency_;
        }

        writeEntry(os, "efficiencies", efficiencies);
    }
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
