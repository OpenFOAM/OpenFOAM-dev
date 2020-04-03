/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::thirdBodyEfficiencies::thirdBodyEfficiencies
(
    const speciesTable& species,
    const scalarList& efficiencies
)
:
    scalarList(efficiencies),
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
    scalarList(species.size()),
    species_(species)
{
    if (dict.found("coeffs"))
    {
        List<Tuple2<word, scalar>> coeffs(dict.lookup("coeffs"));
        if (coeffs.size() != species_.size())
        {
            FatalIOErrorInFunction(dict)
                << "number of efficiencies = " << coeffs.size()
                << " is not equal to the number of species " << species_.size()
                << exit(FatalIOError);
        }

        forAll(coeffs, i)
        {
            operator[](species[coeffs[i].first()]) = coeffs[i].second();
        }
    }
    else
    {
        scalar defaultEff = dict.lookup<scalar>("defaultEfficiency");
        scalarList::operator=(defaultEff);
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::thirdBodyEfficiencies::M(const scalarList& c) const
{
    scalar M = 0;
    forAll(*this, i)
    {
        M += operator[](i)*c[i];
    }

    return M;
}


inline void Foam::thirdBodyEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> coeffs(species_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[i];
        coeffs[i].second() = operator[](i);
    }

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
