/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "dimensionSet.H"
#include "symbols.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dimensionSet::round(const scalar tol)
{
    for (int i=0; i < dimensionSet::nDimensions; ++i)
    {
        scalar integralPart;
        scalar fractionalPart = std::modf(exponents_[i], &integralPart);

        if (mag(fractionalPart - 1.0) <= tol)
        {
            exponents_[i] = 1.0 + integralPart;
        }
        else if (mag(fractionalPart + 1.0) <= tol)
        {
            exponents_[i] = -1.0 + integralPart;
        }
        else if (mag(fractionalPart) <= tol)
        {
            exponents_[i] = integralPart;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dimensionSet::dimensionSet(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::dimensionSet::readNoBegin(Istream& is)
{
    // Read next token
    token nextToken(is);

    if (!nextToken.isNumber())
    {
        // Named dimensions. Parse.
        is.putBack(nextToken);
        reset(symbols::parseNoBegin(is, dimless, dimensions()));
    }
    else
    {
        // Numbered dimensions. Read directly.

        // Read first five dimensions
        exponents_[dimensionSet::MASS] = nextToken.number();
        for (int Dimension=1; Dimension<dimensionSet::CURRENT; Dimension++)
        {
            is >> exponents_[Dimension];
        }

        // Read next token
        token nextToken(is);

        // If next token is another number then read the last two dimensions
        // and then read another token for the end of the dimensionSet
        if (nextToken.isNumber())
        {
            exponents_[dimensionSet::CURRENT] = nextToken.number();
            is >> nextToken;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = nextToken.number();
            is >> nextToken;
        }
        else
        {
            exponents_[dimensionSet::CURRENT] = 0;
            exponents_[dimensionSet::LUMINOUS_INTENSITY] = 0;
        }

        // Check end of dimensionSet
        if (nextToken != token::END_SQR)
        {
            FatalIOErrorInFunction(is)
                << "expected a " << token::END_SQR << " in dimensionSet "
                << endl << "in stream " << is.info()
                << exit(FatalIOError);
        }
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, dimensionSet&)");

    return is;
}


Foam::Istream& Foam::dimensionSet::read(Istream& is)
{
    // Read beginning of dimensionSet
    token startToken(is);

    if (startToken != token::BEGIN_SQR)
    {
        FatalIOErrorInFunction(is)
            << "expected a " << token::BEGIN_SQR << " in dimensionSet"
            << endl << "in stream " << is.info()
            << exit(FatalIOError);
    }

    return readNoBegin(is);
}


Foam::Ostream& Foam::dimensionSet::write(Ostream& os) const
{
    if (dimensionless())
    {
        os << token::BEGIN_SQR << token::END_SQR;
    }
    else
    {
        os << token::BEGIN_SQR;

        for (int d=0; d<dimensionSet::nDimensions-1; d++)
        {
            os << exponents_[d] << token::SPACE;
        }

        os << exponents_[dimensionSet::nDimensions-1] << token::END_SQR;
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const dimensionSet&)");

    return os;
}


void Foam::writeEntry(Ostream& os, const dimensionSet& value)
{
    os << value;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, dimensionSet& dims)
{
    dims.read(is);

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, dimensionSet&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const dimensionSet& dims)
{
    dims.write(os);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const dimensionSet&)");

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<dimensionSet>& ip)
{
    const dimensionSet& dims = ip.t_;

    os << token::BEGIN_SQR;

    for (int first=true, i=0; i<dimensionSet::nDimensions; i++)
    {
        if (mag(dims.exponents_[i]) > dimensionSet::smallExponent)
        {
            if (!first)
            {
                os << token::SPACE;
            }

            os << dimensionSet::dimensionTypeNames_
                  [static_cast<dimensionSet::dimensionType>(i)];

            if (dims.exponents_[i] != 1)
            {
                os << '^' << dims.exponents_[i];
            }

            first = false;
        }
    }

    os << token::END_SQR;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const InfoProxy<dimensionSet>&)");

    return os;
}


// ************************************************************************* //
