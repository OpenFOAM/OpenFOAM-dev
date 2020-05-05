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

#include "thermophysicalFunction.H"
#include "noneThermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalFunction, 0);
    defineRunTimeSelectionTable(thermophysicalFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    const dictionary& dict,
    const word& name
)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing thermophysicalFunction"
            << endl;
    }

    if (dict.isDict(name))
    {
        const dictionary& funcDict(dict.subDict(name));
        const word thermophysicalFunctionType(funcDict.lookup("type"));

        dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(thermophysicalFunctionType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown thermophysicalFunction type "
                << thermophysicalFunctionType
                << nl << nl
                << "Valid thermophysicalFunction types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<thermophysicalFunction>(cstrIter()(funcDict));
    }
    else
    {
        return autoPtr<thermophysicalFunction>
        (
            new thermophysicalFunctions::none(dict.name()/name)
        );
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermophysicalFunction::write(Ostream& os, const word& name) const
{
    os << nl;
    writeKeyword(os, name)
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
    write(os);
    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const thermophysicalFunction& f)
{
    f.write(os);
    return os;
}


// ************************************************************************* //
