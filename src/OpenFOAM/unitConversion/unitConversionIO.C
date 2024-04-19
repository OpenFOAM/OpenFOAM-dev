/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "unitConversion.H"
#include "dictionary.H"
#include "symbols.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unitConversion::unitConversion(Istream& is)
:
    dimensions_(dimless),
    multiplier_(NaN)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::unitConversion::read(const word& keyword, const dictionary& dict)
{
    const unitConversion units(dict.lookup(keyword));

    if (!compare(*this, units, false))
    {
        FatalIOErrorInFunction(dict)
            << "The units " << units.info() << " provided do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    reset(units);
}


bool Foam::unitConversion::readIfPresent
(
    const word& keyword,
    const dictionary& dict
)
{
    const entry* entryPtr = dict.lookupEntryPtr(keyword, false, true);

    if (entryPtr)
    {
        const unitConversion units(entryPtr->stream());

        if (!compare(*this, units, false))
        {
            FatalIOErrorInFunction(dict)
                << "The units " << units.info() << " provided do not match "
                << "the required units " << info()
                << abort(FatalIOError);
        }

        reset(units);

        return true;
    }
    else
    {
        if (dictionary::writeOptionalEntries)
        {
            IOInfoInFunction(dict)
                << "Optional entry '" << keyword << "' is not present,"
                << " the default value '" << info() << "' will be used."
                << endl;
        }

        return false;
    }
}


bool Foam::unitConversion::readIfPresent(Istream& is)
{
    token nextToken(is);
    is.putBack(nextToken);

    if (nextToken != token::BEGIN_SQR) return false;

    unitConversion u(is);

    if (!unitConversion::compare(u, *this, false))
    {
        FatalIOErrorInFunction(is)
            << "The units " << u.info() << " provided do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    if (debug && (any() || !unitConversion::compare(u, *this, true)))
    {
        Info<< "Unit conversion at line " << is.lineNumber()
            << " of file " << is.name()
            << " with factor " << u.multiplier_ << endl;
    }

    reset(u);

    return true;
}


bool Foam::unitConversion::readIfPresent
(
    const word& keyword,
    const dictionary& dict,
    Istream& is
)
{
    token nextToken(is);
    is.putBack(nextToken);

    if (nextToken != token::BEGIN_SQR) return false;

    unitConversion u(is);

    if (!unitConversion::compare(u, *this, false))
    {
        FatalIOErrorInFunction(is)
            << "The units " << u.info() << " of " << keyword
            << " in dictionary " << dict.name() << " do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    if (debug && (any() || !unitConversion::compare(u, *this, true)))
    {
        Info<< "Unit conversion of " << keyword
            << " in dictionary " << dict.name()
            << " with factor " << u.multiplier_ << endl;
    }

    reset(u);

    return true;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, unitConversion& units)
{
    // Read beginning of unitConversion
    token startToken(is);

    if (startToken != token::BEGIN_SQR)
    {
        FatalIOErrorInFunction(is)
            << "expected a " << token::BEGIN_SQR << " in unitConversion"
            << endl << "in stream " << is.info()
            << exit(FatalIOError);
    }

    // Peek at the next token
    token nextToken(is);
    is.putBack(nextToken);

    if (!nextToken.isNumber())
    {
        // Named units. Parse.
        units.reset(symbols::parseNoBegin(is, unitless, Foam::units()));
    }
    else
    {
        // Read the dimensions
        units.dimensions_.readNoBegin(is);

        // Read the dimensionless units if present, or set to zero
        token nextToken;
        if (!is.eof())
        {
            is  >> nextToken;
        }
        if (nextToken == token::BEGIN_SQR)
        {
            for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
            {
                is  >> units.exponents_[i];
            }

            // Check end of dimensionless units
            token endToken(is);
            if (endToken != token::END_SQR)
            {
                FatalIOErrorInFunction(is)
                    << "expected a " << token::END_SQR
                    << " in unitConversion " << endl << "in stream "
                    << is.info() << exit(FatalIOError);
            }

            // Read the multiplier if present, or set to unity
            token nextToken;
            if (!is.eof())
            {
                is  >> nextToken;
            }
            if (nextToken == token::BEGIN_SQR)
            {
                is  >> units.multiplier_;

                // Check end of multiplier
                token endToken(is);
                if (endToken != token::END_SQR)
                {
                    FatalIOErrorInFunction(is)
                        << "expected a " << token::END_SQR
                        << " in unitConversion " << endl << "in stream "
                        << is.info() << exit(FatalIOError);
                }
            }
            else
            {
                units.multiplier_ = 1;
                if (!nextToken.undefined())
                {
                    is.putBack(nextToken);
                }
            }
        }
        else
        {
            for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
            {
                units.exponents_[i] = 0;
            }
            units.multiplier_ = 1;
            if (!nextToken.undefined())
            {
                is.putBack(nextToken);
            }
        }
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, unitConversion&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const unitConversion& units)
{
    // Write the dimensions
    os  << units.dimensions_;

    // Determine if any dimensionless units are non-zero
    bool nonZeroDimlessUnits = false;
    for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
    {
        nonZeroDimlessUnits =
            nonZeroDimlessUnits
         || mag(units.exponents_[i]) > unitConversion::smallExponent;
    }

    // Determine if the multiplier is non-unity
    bool nonUnityMultiplier = units.multiplier_ != 1;

    // Write the dimensionless units if any are non-zero or we have a
    // multiplier to write out afterwards
    if (nonZeroDimlessUnits || nonUnityMultiplier)
    {
        os  << token::BEGIN_SQR;
        for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
        {
            if (i) os  << token::SPACE;
            os  << units.exponents_[i];
        }
        os  << token::END_SQR;
    }

    // Write the multiplier if it is non-unity
    if (nonUnityMultiplier)
    {
        os  << token::BEGIN_SQR << units.multiplier_ << token::END_SQR;
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const unitConversion&)");

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<unitConversion>& ip
)
{
    const unitConversion& units = ip.t_;

    // Filter out special cases
    if (units.any())
    {
        return os << token::BEGIN_SQR << "<any>" << token::END_SQR;
    }
    if (units.none())
    {
        return os << token::BEGIN_SQR << "<none>" << token::END_SQR;
    }

    // Write the dimensions
    os << units.dimensions_.info();

    // Determine if any dimensionless units are non-zero
    bool nonZeroDimlessUnits = false;
    for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
    {
        nonZeroDimlessUnits =
            nonZeroDimlessUnits
         || mag(units.exponents_[i]) > unitConversion::smallExponent;
    }

    // Determine if the multiplier is non-unity
    bool nonUnityMultiplier = units.multiplier_ != 1;

    // Write the dimensionless units if any are non-zero or we have a
    // multiplier to write out afterwards
    if (nonZeroDimlessUnits || nonUnityMultiplier)
    {
        // Write the dimensionless units
        os << token::BEGIN_SQR;

        for (int first=true, i=0; i<unitConversion::nDimlessUnits; i++)
        {
            if (mag(units.exponents_[i]) > unitConversion::smallExponent)
            {
                if (!first)
                {
                    os << token::SPACE;
                }

                os << unitConversion::dimlessUnitTypeNames_
                      [static_cast<unitConversion::dimlessUnitType>(i)];

                if (units.exponents_[i] != 1)
                {
                    os << '^' << units.exponents_[i];
                }

                first = false;
            }
        }

        os << token::END_SQR;
    }

    // Write the multiplier if it is non-unity
    if (nonUnityMultiplier)
    {
        os  << token::BEGIN_SQR << units.multiplier_ << token::END_SQR;
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const InfoProxy<unitConversion>&)");

    return os;
}


// ************************************************************************* //
