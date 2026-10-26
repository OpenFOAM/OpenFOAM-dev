/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "reactionCoeffs.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reactionCoeffs::setSpecieIndices()
{
    // Fill specie indices with left then right hand side indices
    specieIndices_.resize(lhs_.size() + rhs_.size());
    forAll(lhs_, i)
    {
        specieIndices_[i] = lhs_[i].index;
    }
    forAll(rhs_, i)
    {
        specieIndices_[i + lhs_.size()] = rhs_[i].index;
    }

    // Sort
    Foam::sort(specieIndices_);

    // Make unique
    label i0 = 0;
    for (label i = 1; i < specieIndices_.size(); ++ i)
    {
        if (specieIndices_[i] != specieIndices_[i0]) i0 ++;
        specieIndices_[i0] = specieIndices_[i];
    }
    specieIndices_.resize(i0 + 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionCoeffs::reactionCoeffs
(
    const word& name,
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    name_(name),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
{
    setSpecieIndices();
}


Foam::reactionCoeffs::reactionCoeffs
(
    const speciesTable& species,
    const dictionary& dict
)
:
    name_(dict.dictName()),
    species_(species)
{
    // Read as a string
    ITstream& its = dict.lookup("reaction");
    const string str(its);

    // Parse the string
    IStringStream iss(str);
    bool right = false;
    DynamicList<specieCoeffs> dlrhs;
    for (;;)
    {
        dlrhs.append(specieCoeffs(species_, iss));

        if (!iss.good()) break;

        token t(iss);

        if (t == token::ADD) continue;

        if (t == token::ASSIGN)
        {
            if (!right)
            {
                lhs_.transfer(dlrhs);
                right = true;
                continue;
            }
            else
            {
                FatalIOErrorInFunction(its)
                    << "Reaction has multiple '=' operators"
                    << exit(FatalIOError);
            }
        }

        FatalIOErrorInFunction(its)
            << "Reaction format not recognised"
            << exit(FatalIOError);
    }

    if (!right)
    {
        FatalIOErrorInFunction(its)
            << "Reaction has no '=' operator"
            << exit(FatalIOError);
    }

    rhs_.transfer(dlrhs);

    setSpecieIndices();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reactionCoeffs::write(Ostream& os) const
{
    OStringStream oss;

    forAll(lhs_, i)
    {
        if (i > 0) oss << " + ";
        lhs_[i].write(species_, oss);
    }

    oss << " = ";

    forAll(rhs_, i)
    {
        if (i > 0) oss << " + ";
        rhs_[i].write(species_, oss);
    }

    writeEntry(os, "reaction", oss.str());
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const reactionCoeffs& r)
{
    r.write(os);
    return os;
}


// ************************************************************************* //
