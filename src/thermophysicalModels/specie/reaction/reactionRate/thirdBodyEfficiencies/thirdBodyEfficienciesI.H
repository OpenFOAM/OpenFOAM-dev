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
    const scalarList& efficiencies
)
:
    scalarField(efficiencies),
    species_(species)
{
    if (size() != species_.size())
    {
        FatalErrorInFunction
            << "number of efficiencies = " << size()
            << " is not equal to the number of species " << species_.size()
            << exit(FatalError);
    }
}


inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const dictionary& dict
)
:
    scalarField(species.size(), 1.0),
    species_(species)
{
    if (dict.found("coeffs"))
    {
        List<Tuple2<word, scalar>> coeffs(dict.lookup("coeffs"));
        forAll(coeffs, i)
        {
            operator[](species[coeffs[i].first()]) = coeffs[i].second();
        }
    }
    else
    {
        scalar defaultEff = dict.lookup<scalar>("defaultEfficiency");
        scalarField::operator=(defaultEff);
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M
(
    const reaction::state& s
) const
{
    scalar M = 0;
    forAll(*this, i)
    {
        M += operator[](i)*s.c[i];
    }

    return M;
}


inline const Foam::scalarField&
Foam::thirdBodyEfficiencies::dMdc() const
{
    return *this;
}


inline void Foam::thirdBodyEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> coeffs(species_.size());
    label n = 0;
    forAll(coeffs, i)
    {
        if (operator[](i) != 1)
        {
            coeffs[n].first() = species_[i];
            coeffs[n].second() = operator[](i);
            n++;
        }
    }
    coeffs.setSize(n);

    writeEntry(os, "coeffs", coeffs);
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
