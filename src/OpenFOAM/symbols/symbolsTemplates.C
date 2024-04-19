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

#include "symbols.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::symbols::parseNoBegin
(
    const label lastPrior,
    tokeniser& tis,
    const Type& identity,
    const HashTable<Type>& table
)
{
    Type result(identity);

    // Get initial token
    token nextToken(tis.nextToken());

    // Store type of last token read. Used to detect two consecutive
    // symbols and assume multiplication.
    bool haveReadSymbol = false;

    while (true)
    {
        if (nextToken.isWord())
        {
            // Named unit conversion. Multiply.
            result.reset(result*table[nextToken.wordToken()]);
            haveReadSymbol = true;
        }
        else if (nextToken.isNumber())
        {
            // A number. This makes no sense.
            FatalIOErrorInFunction(tis.stream())
                << "Illegal token " << nextToken << exit(FatalIOError);
        }
        else if (nextToken.isPunctuation())
        {
            const label nextPrior = tokeniser::priority(nextToken);

            if (nextToken.pToken() == token::BEGIN_SQR)
            {
                // Start another set of symbols? This makes no sense.
                FatalIOErrorInFunction(tis.stream())
                    << "Illegal token " << nextToken << exit(FatalIOError);
                return result;
            }
            else if (nextToken.pToken() == token::END_SQR)
            {
                // End the units
                tis.putBack(nextToken);
                return result;
            }
            else if (nextToken.pToken() == token::BEGIN_LIST)
            {
                // Parenthesis. Evaluate the sub-units and multiply.
                result.reset
                (
                    result*parseNoBegin(nextPrior, tis, identity, table)
                );

                // Check that the parentheses end
                token t = tis.nextToken();
                if (!t.isPunctuation() || t.pToken() != token::END_LIST)
                {
                    FatalIOErrorInFunction(tis.stream())
                        << "Illegal token " << t << exit(FatalIOError);
                }

                haveReadSymbol = true;
            }
            else if (nextToken.pToken() == token::END_LIST)
            {
                // End the sub-units
                tis.putBack(nextToken);
                return result;
            }
            else if (nextToken.pToken() == token::MULTIPLY)
            {
                // Multiply operator
                if (nextPrior > lastPrior)
                {
                    // This has priority. Evaluate the next units and multiply.
                    result.reset
                    (
                        result*parseNoBegin(nextPrior, tis, identity, table)
                    );
                }
                else
                {
                    // Restore the token
                    tis.putBack(nextToken);
                    return result;
                }

                haveReadSymbol = false;
            }
            else if (nextToken.pToken() == token::DIVIDE)
            {
                // Divide operator. As above.
                if (nextPrior > lastPrior)
                {
                    result.reset
                    (
                        result/parseNoBegin(nextPrior, tis, identity, table)
                    );
                }
                else
                {
                    tis.putBack(nextToken);
                    return result;
                }

                haveReadSymbol = false;
            }
            else if (nextToken.pToken() == '^')
            {
                // Power operator
                if (nextPrior > lastPrior)
                {
                    token t = tis.nextToken();
                    if (!t.isScalar())
                    {
                        FatalIOErrorInFunction(tis.stream())
                            << "Invalid power " << t << exit(FatalIOError);
                    }
                    result.reset(pow(result, t.scalarToken()));
                }
                else
                {
                    tis.putBack(nextToken);
                    return result;
                }

                haveReadSymbol = true;
            }
            else
            {
                FatalIOErrorInFunction(tis.stream())
                    << "Illegal token " << nextToken << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction(tis.stream())
                << "Illegal token " << nextToken << exit(FatalIOError);
        }

        if (!tis.hasToken())
        {
            break;
        }

        nextToken = tis.nextToken();
        if (nextToken.error())
        {
            break;
        }

        if (haveReadSymbol && (nextToken.isWord() || nextToken.isNumber()))
        {
            // Consecutive symbols are multiplied
            tis.putBack(nextToken);
            nextToken = token(token::MULTIPLY);
        }
    }

    return result;
}


template<class Type>
Type Foam::symbols::parseNoBegin
(
    Istream& is,
    const Type& identity,
    const HashTable<Type>& table
)
{
    symbols::tokeniser tis(is);

    return parseNoBegin(0, tis, identity, table);
}


// ************************************************************************* //
