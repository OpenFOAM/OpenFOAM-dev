/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
            << "The units " << units.info() << " of " << keyword
            << " in dictionary " << dict.name() << " do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    reset(units);
}


void Foam::unitConversion::read(Istream& is)
{
    const unitConversion units(is);

    if (!compare(*this, units, false))
    {
        FatalIOErrorInFunction(is)
            << "The units " << units.info() << " provided do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    reset(units);
}


void Foam::unitConversion::read
(
    const word& keyword,
    const dictionary& dict,
    Istream& is
)
{
    const unitConversion units(is);

    if (!compare(*this, units, false))
    {
        FatalIOErrorInFunction(dict)
            << "The units " << units.info() << " of " << keyword
            << " in dictionary " << dict.name() << " do not match "
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
                << "The units " << units.info() << " of " << keyword
                << " in dictionary " << dict.name() << " do not match "
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

    const unitConversion units(is);

    if (!unitConversion::compare(units, *this, false))
    {
        FatalIOErrorInFunction(is)
            << "The units " << units.info() << " provided do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    if (debug && (any() || !unitConversion::compare(units, *this, true)))
    {
        Info<< "Unit conversion at line " << is.lineNumber()
            << " of file " << is.name()
            << " with factor " << units.multiplier_ << endl;
    }

    reset(units);

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

    const unitConversion units(is);

    if (!unitConversion::compare(units, *this, false))
    {
        FatalIOErrorInFunction(dict)
            << "The units " << units.info() << " of " << keyword
            << " in dictionary " << dict.name() << " do not match "
            << "the required units " << info()
            << abort(FatalIOError);
    }

    if (debug && (any() || !unitConversion::compare(units, *this, true)))
    {
        Info<< "Unit conversion of " << keyword
            << " in dictionary " << dict.name()
            << " with factor " << units.multiplier_ << endl;
    }

    reset(units);

    return true;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, unitConversion& units)
{
    token nextToken;

    // Read the next delimiting token. This must be the start bracket.
    is >> nextToken;
    if (nextToken != token::BEGIN_SQR)
    {
        FatalIOErrorInFunction(is)
            << "expected a " << token::BEGIN_SQR << " in unitConversion"
            << endl << "in stream " << is.info() << ", got a "
            << nextToken << exit(FatalIOError);
    }

    // Peek at the next token
    is >> nextToken;
    is.putBack(nextToken);

    // If not a number or separator, then these are named units. Parse.
    if (!nextToken.isNumber() && nextToken != token::COLON)
    {
        // Named units. Parse. Note: Use an explicit construction of the
        // identity unit instead of 'unitless' because this function may be
        // used in a static context before 'unitless' is available.
        units.reset
        (
            symbols::parseNoBeginOrEnd
            (
                is,
                unitConversion(dimless, 0, 0, 1),
                Foam::units()
            )
        );

        // Read the next delimiting token. This must be the end bracket.
        is >> nextToken;
        if (nextToken != token::END_SQR)
        {
            FatalIOErrorInFunction(is)
                << "expected a " << token::END_SQR << " in unitConversion "
                << endl << "in stream " << is.info() << ", got a "
                << nextToken << exit(FatalIOError);
        }

        // Check state of Istream
        is.check("Istream& operator>>(Istream&, unitConversion&)");

        return is;
    }

    // Otherwise these are numbered units. Read directly...

    // Read the dimensions
    units.dimensions_.readNoBeginOrEnd(is);

    // Read the next delimiting token. If a separator, then there are
    // dimensionless units and a multiplier to read. If it is an end bracket,
    // then the dimensionless units are zero and the multiplier is one and the
    // parsing is finished. Otherwise the parsing has failed.
    is >> nextToken;
    if (nextToken == token::COLON)
    {
        // Peek at the next token
        is >> nextToken;
        is.putBack(nextToken);

        // Read the dimensionless units if present, or set to zero
        if (!nextToken.isNumber())
        {
            for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
            {
                units.exponents_[i] = 0;
            }
        }
        else
        {
            for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
            {
                is  >> units.exponents_[i];
            }
        }

        // Read the next delimiting token. If a separator then there is a
        // multiplier to read. If it is an end bracket then the multiplier is
        // one and the parsing is finished. Otherwise the parsing has failed.
        is >> nextToken;
        if (nextToken == token::COLON)
        {
            // Peek at the next token
            is >> nextToken;
            is.putBack(nextToken);

            // Read the multiplier if present, or set to unity
            if (!nextToken.isNumber())
            {
                units.multiplier_ = 1;
            }
            else
            {
                is >> units.multiplier_;
            }

            // Read the next delimiting token. This must be the end bracket.
            is >> nextToken;
            if (nextToken != token::END_SQR)
            {
                FatalIOErrorInFunction(is)
                    << "expected a " << token::END_SQR << " in unitConversion "
                    << endl << "in stream " << is.info() << ", got a "
                    << nextToken << exit(FatalIOError);
            }
        }
        else if (nextToken == token::END_SQR)
        {
            units.multiplier_ = 1;
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "expected a " << token::END_SQR << " or a " << token::COLON
                << " in unitConversion " << endl << "in stream " << is.info()
                << ", got a " << nextToken << exit(FatalIOError);
        }
    }
    else if (nextToken == token::END_SQR)
    {
        for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
        {
            units.exponents_[i] = 0;
        }

        units.multiplier_ = 1;
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "expected a " << token::END_SQR << " or a " << token::COLON
            << " in unitConversion " << endl << "in stream " << is.info()
            << ", got a " << nextToken << exit(FatalIOError);
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, unitConversion&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const unitConversion& units)
{
    // Write the start
    os << token::BEGIN_SQR;

    // Write the dimensions
    units.dimensions_.writeNoBeginOrEnd(os);

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

    // Write a separator if there is anything to follow
    if (nonZeroDimlessUnits || nonUnityMultiplier)
    {
        os  << token::SPACE << token::COLON;
    }

    // Write the dimensionless units if any are non-zero
    if (nonZeroDimlessUnits)
    {
        for (int i = 0; i < unitConversion::nDimlessUnits; ++ i)
        {
            os  << token::SPACE << units.exponents_[i];
        }
    }

    // Write a separator if there is anything to follow
    if (nonUnityMultiplier)
    {
        os  << token::SPACE << token::COLON;
    }

    // Write the multiplier if it is non-unity
    if (nonUnityMultiplier)
    {
        os  << token::SPACE << units.multiplier_;
    }

    // Write the end
    os  << token::END_SQR;

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

    // Write the start
    os << token::BEGIN_SQR;

    // Write the dimensions
    units.dimensions_.writeInfoNoBeginOrEnd(os);

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

    // Write a separator if there is anything to follow
    if (nonZeroDimlessUnits || nonUnityMultiplier)
    {
        os  << token::SPACE << token::COLON;
    }

    // Write the dimensionless units if any are non-zero
    if (nonZeroDimlessUnits)
    {
        for (int i=0; i<unitConversion::nDimlessUnits; i++)
        {
            if (mag(units.exponents_[i]) > unitConversion::smallExponent)
            {
                os << token::SPACE << unitConversion::dimlessUnitTypeNames_
                      [static_cast<unitConversion::dimlessUnitType>(i)];

                if (units.exponents_[i] != 1)
                {
                    os << '^' << units.exponents_[i];
                }
            }
        }
    }

    // Write a separator if there is anything to follow
    if (nonUnityMultiplier)
    {
        os  << token::SPACE << token::COLON;
    }

    // Write the multiplier if it is non-unity
    if (nonUnityMultiplier)
    {
        os  << token::SPACE << units.multiplier_;
    }

    // Write the end
    os  << token::END_SQR;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const InfoProxy<unitConversion>&)");

    return os;
}


// ************************************************************************* //
