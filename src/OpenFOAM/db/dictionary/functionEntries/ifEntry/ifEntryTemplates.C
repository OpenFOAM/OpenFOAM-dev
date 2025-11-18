/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "Switch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Context>
bool Foam::functionEntries::ifEntry::execute
(
    DynamicList<filePos>& stack,
    const dictionary& contextDict,
    Context& context,
    Istream& is
) const
{
    const label nNested = stack.size();

    stack.append(filePos(is.name(), is.lineNumber()));

    const string arg
    (
        operator[](1).stringToken() + char(token::END_STATEMENT)
    );
    IStringStream argStream(arg);
    argStream.lineNumber() = operator[](1).lineNumber();
    const primitiveEntry e("ifEntry", contextDict, argStream);
    const Switch doIf(e.stream());

    bool ok = ifeqEntry::execute(doIf, stack, contextDict, context, is);

    if (stack.size() != nNested)
    {
        FatalIOErrorInFunction(contextDict)
            << "Did not find matching #endif for "
            << typeName << " condition starting"
            << " at line " << stack.last().second()
            << " in file " <<  stack.last().first() << exit(FatalIOError);
    }

    return ok;
}


// ************************************************************************* //
