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
    fvSolutionCombine

Description
    Simple utility for combining fvSolution solution entries.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "wordRe.H"
#include "OSspecific.H"

using namespace Foam;

// check for identical dictionary content, regardless of the order
bool checkDictionaryContent(const dictionary& dict1, const dictionary& dict2)
{
    // trivial cases first
    if (&dict1 == &dict2)
    {
        return true;
    }
    else if (dict1.size() != dict2.size())
    {
        return false;
    }


    forAllConstIter(dictionary, dict1, iter1)
    {
        const entry* entryPtr = dict2.lookupEntryPtr
        (
            iter1().keyword(),
            false,
            false
        );

        if (!entryPtr)
        {
            return false;
        }

        const entry& entry1 = iter1();
        const entry& entry2 = *entryPtr;

        bool ok = false;
        if (entry1.isDict())
        {
            if (entry2.isDict())
            {
                ok = checkDictionaryContent(entry1.dict(), entry2.dict());
            }
        }
        else
        {
            ok = (entry1 == entry2);
        }


        if (!ok)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption("rewrite");
    argList::addBoolOption("show");

    argList args(argc, argv);

    Time runTime(args.rootPath(), args.caseName());

    const word dictName("fvSolution");

    bool optRewrite = args.optionFound("rewrite");
    bool optShow = args.optionFound("show");

    IOdictionary solutionDict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!solutionDict.found("solvers"))
    {
        Info<<"no solvers entry found in : " << dictName << endl;
        return 2;
    }

    if (optRewrite && solutionDict.instance() != runTime.system())
    {
        Info<<"instance is not " << runTime.system()
            << "- disabling rewrite for this file" << nl;
        optRewrite = false;
    }

    dictionary& solverDict = solutionDict.subDict("solvers");

    wordList names = solverDict.toc();
    wordList oldNames = names;

    bool changed = false;
    for (label orig = 0; orig < names.size()-1; ++orig)
    {
        // skip patterns or entries that have already been done
        if (names[orig].empty() || wordRe::isPattern(names[orig]))
        {
            continue;
        }

        const dictionary& dict1 = solverDict.subDict(names[orig]);

        for (label check = orig+1; check < names.size(); ++check)
        {
            // skip patterns or entries that have already been done
            if (names[check].empty() || wordRe::isPattern(names[check]))
            {
                continue;
            }

            const dictionary& dict2 = solverDict.subDict(names[check]);

            // check for identical content
            if (checkDictionaryContent(dict1, dict2))
            {
                names[orig] += "|" + names[check];
                names[check].clear();
                changed = true;
            }
        }
    }

    if (changed)
    {
        forAll(names, nameI)
        {
            if (names[nameI].empty())
            {
                solverDict.remove(oldNames[nameI]);
                Info<<"  #remove " << oldNames[nameI];
            }
            else
            {
                Info<< "  " << oldNames[nameI];

                if (names[nameI] != oldNames[nameI])
                {
                    // make "(abc|def)" pattern
                    keyType renamed( "(" + names[nameI] + ")", true);

                    solverDict.changeKeyword(oldNames[nameI], renamed);

                    Info<< " -> " << renamed;
                }
            }
            Info<< endl;
        }

        if (optRewrite)
        {
            mvBak(solutionDict.objectPath(), "orig");
            Info<< "Backup to .orig" << nl
                << "Writing " << solutionDict.objectPath() << nl << endl;

            solutionDict.regIOobject::write();
        }
        else if (optShow)
        {
            IOobject::writeDivider(Info);
            solutionDict.dictionary::write(Info, false);
        }
        else
        {
            Info<< "\nFile not rewritten" << endl;
        }
    }
    else
    {
        Info<< "no changes" << endl;
    }

    Info<< "Done\n" << endl;

    return changed ? 0 : 1;
}


// ************************************************************************* //
