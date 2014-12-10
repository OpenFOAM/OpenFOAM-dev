/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "error.H"

#include "thermophysicalFunction.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalFunction, 0);
    defineRunTimeSelectionTable(thermophysicalFunction, Istream);
    defineRunTimeSelectionTable(thermophysicalFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    Istream& is
)
{
    if (debug)
    {
        Info<< "thermophysicalFunction::New(Istream&) : "
            << "constructing thermophysicalFunction"
            << endl;
    }

    const word thermophysicalFunctionType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(thermophysicalFunctionType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("thermophysicalFunction::New(Istream&)")
            << "Unknown thermophysicalFunction type "
            << thermophysicalFunctionType
            << nl << nl
            << "Valid thermophysicalFunction types are :" << endl
            << IstreamConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<thermophysicalFunction>(cstrIter()(is));
}


Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "thermophysicalFunction::New(const dictionary&) : "
            << "constructing thermophysicalFunction"
            << endl;
    }

    const word thermophysicalFunctionType(dict.lookup("functionType"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermophysicalFunctionType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("thermophysicalFunction::New(const dictionary&)")
            << "Unknown thermophysicalFunction type "
            << thermophysicalFunctionType
            << nl << nl
            << "Valid thermophysicalFunction types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<thermophysicalFunction>(cstrIter()(dict));
}


// ************************************************************************* //
