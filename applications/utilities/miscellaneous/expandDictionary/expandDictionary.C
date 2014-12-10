/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    expandDictionary

Description
    Read the dictionary provided as an argument, expand the macros etc. and
    write the resulting dictionary to standard output.

Usage
    - expandDictionary inputDict [OPTION]

    \param -list \n
    Report the #include/#includeIfPresent to stdout only.

Note
    The \c -list option can be useful when determining which files
    are actually included by a directory. It can also be used to
    determine which files may need to be copied when transferring
    simulation to another environment. The following code snippet
    could be a useful basis for such cases:

    \verbatim
        for i in . 0 constant system
        do
            find $i -maxdepth 1 -type f -exec expandDictionary -list '{}' \;
        done | sed -ne '/^"\//!{ s/^"//; s/"$//; p }' | sort | uniq
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "IOobject.H"
#include "dictionary.H"
#include "includeEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read the specified dictionary file, expand the macros etc. and write\n"
        "the resulting dictionary to standard output."
    );

    argList::addBoolOption
    (
        "list",
        "Report the #include/#includeIfPresent to stdout only"
    );

    argList::noBanner();
    argList::noParallel();
    argList::validArgs.append("inputDict");
    argList args(argc, argv);

    const string dictName = args[1];

    const bool listOpt = args.optionFound("list");

    if (listOpt)
    {
        Foam::functionEntries::includeEntry::report = true;
    }

    dictionary dict(IFstream(dictName)(), true);

    if (!listOpt)
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictName << "\n//\n";
        dict.write(Info, false);
        IOobject::writeDivider(Info);
    }

    return 0;
}


// ************************************************************************* //
