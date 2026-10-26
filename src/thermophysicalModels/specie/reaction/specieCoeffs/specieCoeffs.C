/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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

#include "specieCoeffs.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specieCoeffs::specieCoeffs()
:
    index(-labelMax),
    stoichCoeff(NaN),
    exponent()
{}


Foam::specieCoeffs::specieCoeffs(const speciesTable& species, Istream& is)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
    }
    else
    {
        stoichCoeff = 1;
        is.putBack(t);
    }

    exponent = stoichCoeff;

    word specieName(is);

    const size_t i = specieName.find('^');

    if (i != word::npos)
    {
        const string exponentStr =
            specieName(i + 1, specieName.size() - i - 1);

        exponent = scalar(atof(exponentStr.c_str()));

        if (scalar(exponent) < 0)
        {
            FatalIOErrorInFunction(is)
                << "Negative exponent for specie " << specieName
                << " is not supported."
                << exit(FatalIOError);
        }

        specieName = specieName(0, i);
    }

    if (!species.found(specieName))
    {
        FatalIOErrorInFunction(is)
            << "Specie " << specieName
            << " not found in table " << species
            << exit(FatalIOError);
    }

    index = species[specieName];
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::specieCoeffs::write(const speciesTable& species, Ostream& os) const
{
    if (mag(stoichCoeff - 1) > small)
    {
        os << stoichCoeff;
    }

    os << species[index];

    if (mag(exponent.operator scalar() - stoichCoeff) > small)
    {
        os << '^' << exponent;
    }
}


// ************************************************************************* //
