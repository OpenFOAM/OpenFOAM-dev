/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "DimensionedFunction1.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Dimensioned<Type>::read
(
    const word& name,
    const dictionary& dict
)
{
    // Read units from the given stream
    auto readUnits = []
    (
        Istream& is,
        const word& prefix,
        dimensionedScalar& units
    )
    {
        units.name() =
            prefix + (prefix.empty() ? 'u' : 'U') + "nits";

        dimensionSet d(dimless);
        d.read(is, units.value());

        if (d != units.dimensions())
        {
            FatalIOErrorInFunction(is)
                << "The " << prefix << (prefix.empty() ? "" : "-")
                << "dimensions " << d << " provided do not match the required "
                << "dimensions " << units.dimensions()
                << exit(FatalIOError);
        }
    };

    // Read units if present in the given dictionary
    auto lookupUnitsIfPresent = [&readUnits]
    (
        const dictionary& dict,
        const word& prefix,
        dimensionedScalar& units
    )
    {
        const word unitsKey =
            prefix + (prefix.empty() ? 'u' : 'U') + "nits";
        const word dimsKey =
            prefix + (prefix.empty() ? 'd' : 'D') + "imensions";

        const bool haveUnits = dict.found(unitsKey);
        const bool haveDims = dict.found(dimsKey);

        if (haveUnits && haveDims)
        {
            FatalIOErrorInFunction(dict)
                << "Both " << unitsKey << " and " << dimsKey
                << " are specified. Only one is permitted."
                << exit(FatalError);
        }

        if (haveUnits)
        {
            units = dimensionedScalar(unitsKey, units.dimensions(), dict);
        }

        if (haveDims)
        {
            readUnits(dict.lookup(dimsKey), prefix, units);
        }
    };

    // If the function is a dictionary (preferred) then read straightforwardly
    if (dict.isDict(name))
    {
        const dictionary& coeffsDict(dict.subDict(name));

        lookupUnitsIfPresent(coeffsDict, "x", xUnits_);
        lookupUnitsIfPresent(coeffsDict, "", units_);

        value_.reset(Function1<Type>::New(name, dict).ptr());

        return;
    }

    // Find the entry
    Istream& is(dict.lookup(name));

    // Peek at the first token
    token firstToken(is);
    is.putBack(firstToken);

    // Read the type, or assume constant
    const word Function1Type =
        firstToken.isWord() ? word(is) : Constant<Type>::typeName;

    // If the entry is not a type followed by a end statement then
    // construct the function from the stream
    if (!firstToken.isWord() || !is.eof())
    {
        // Peek at the next token
        token nextToken(is);
        is.putBack(nextToken);

        // Read dimensions if they are provided
        if (nextToken == token::BEGIN_SQR)
        {
            readUnits(is, "", units_);
        }

        // Construct from a stream
        value_.reset(Function1<Type>::New(name, Function1Type, is).ptr());

        return;
    }

    // Otherwise, construct from the current or coeffs dictionary
    const dictionary& coeffsDict =
        dict.optionalSubDict(name + "Coeffs");

    lookupUnitsIfPresent(coeffsDict, "x", xUnits_);
    lookupUnitsIfPresent(coeffsDict, "", units_);

    value_.reset(Function1<Type>::New(name, dict).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Dimensioned<Type>::Dimensioned
(
    const word& name,
    const dimensionSet& xDimensions,
    const dimensionSet& dimensions,
    const dictionary& dict
)
:
    FieldFunction1<Type, Dimensioned<Type>>(name),
    xUnits_(word::null, xDimensions, scalar(1)),
    units_(word::null, dimensions, scalar(1)),
    value_(nullptr)
{
    read(name, dict);
}


template<class Type>
Foam::Function1s::Dimensioned<Type>::Dimensioned(const Dimensioned<Type>& df1)
:
    FieldFunction1<Type, Dimensioned<Type>>(df1),
    xUnits_(df1.xUnits_),
    units_(df1.units_),
    value_(df1.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Dimensioned<Type>::~Dimensioned()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Dimensioned<Type>::write(Ostream& os) const
{
    if (xUnits_.name() != word::null)
    {
        writeKeyword(os, xUnits_.name());
        writeEntry(os, xUnits_.dimensions());
        os << token::SPACE;
        writeEntry(os, xUnits_.value());
        os << token::END_STATEMENT << endl;
    }
    if (units_.name() != word::null)
    {
        writeKeyword(os, units_.name());
        writeEntry(os, units_.dimensions());
        os << token::SPACE;
        writeEntry(os, units_.value());
        os << token::END_STATEMENT << endl;
    }
    value_().write(os);
}


// ************************************************************************* //
