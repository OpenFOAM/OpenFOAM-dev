/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "helpType.H"
#include "doxygenXmlParser.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helpType, 0);
    defineRunTimeSelectionTable(helpType, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::helpType::doxygenPath() const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));

    label dirI = -1;
    forAll(docDirs, i)
    {
        if (isDir(docDirs[i].expand()))
        {
            dirI = i;
            break;
        }
    }

    if (dirI == -1)
    {
        Info<< "No Doxygen sources found under search paths: "
            << docDirs << endl;
        return fileName();
    }

    return docDirs[dirI];
}


void Foam::helpType::displayDocOptions
(
    const string& searchStr,
    const bool exactMatch
) const
{
    fileName doxyPath(doxygenPath());

    if (doxyPath.empty())
    {
        return;
    }

    Info<< "Found doxygen help in " << doxyPath.c_str() << nl << endl;

    doxygenXmlParser parser
    (
        doxyPath/"../DTAGS",
        "tagfile",
        searchStr,
        exactMatch
    );

    Info<< "Valid types include:" << nl << SortableList<word>(parser.toc());
}


void Foam::helpType::displayDoc
(
    const word& className,
    const string& searchStr,
    const bool exactMatch
) const
{
    fileName doxyPath(doxygenPath());

    if (doxyPath.empty())
    {
        return;
    }

    Info<< "Found doxygen help in " << doxyPath.c_str() << nl << endl;

    string docBrowser = getEnv("FOAM_DOC_BROWSER");
    if (docBrowser.empty())
    {
        const dictionary& docDict =
            debug::controlDict().subDict("Documentation");
        docDict.lookup("docBrowser") >> docBrowser;
    }

    doxygenXmlParser parser
    (
        doxyPath/"../DTAGS",
        "tagfile",
        searchStr,
        exactMatch
    );

    if (debug)
    {
        Info<< parser;
    }

    if (parser.found(className))
    {
        fileName docFile(doxyPath/parser.subDict(className).lookup("filename"));

        // can use FOAM_DOC_BROWSER='application file://%f' if required
        docBrowser.replaceAll("%f", docFile);

        fileName classDirectory(parser.subDict(className).lookup("path"));
        word classFile(parser.subDict(className).lookup("name"));

        Info<< "Showing documentation for type " << className << nl << endl;

        Info<< "Source file: " << classDirectory.c_str() << classFile << nl
            << endl;

        system(docBrowser);
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::helpType::displayDoc"
            "("
                "const word&, "
                "const string&, "
                "const bool"
            ")"
        )
            << "No help for type " << className << " found."
            << "  Valid options include:" << SortableList<word>(parser.toc())
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpType::helpType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpType::~helpType()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpType::init()
{
    argList::addOption
    (
        "browse",
        "word",
        "display documentation for boundary condition in browser"
    );
}


// ************************************************************************* //
