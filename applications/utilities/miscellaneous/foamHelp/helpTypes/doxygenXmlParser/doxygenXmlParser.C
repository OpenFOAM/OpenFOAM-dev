/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "doxygenXmlParser.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::doxygenXmlParser::doxygenXmlParser
(
    const fileName& fName,
    const string& startTag,
    const string& searchStr,
    const bool exactMatch
)
:
    dictionary(dictionary::null)
{
    // Pre-construct and compile regular expressions
    const wordRe nameRe(".*.H", wordRe::DETECT);
    const wordRe searchStrRe(searchStr, wordRe::DETECT);

    // Pre-construct constant strings and names to speed-up comparisons
    const string slashStartTag('/' + startTag);
    const string kindFileStr("kind=\"file\"");
    const word compoundWord("compound");
    const word nameWord("name");
    const word pathWord("path");
    const word filenameWord("filename");

    IFstream is(fName);

    // Skip forward to entry name
    skipForward(is, startTag);

    char c;

    while (is.get(c))
    {
        if (c == '<')
        {
            // If in block, read block name
            string blockName;
            string params;
            bool readingParam = false;
            while (is.get(c)  && c != '>')
            {
                if (c == ' ')
                {
                    readingParam = true;
                }
                else
                {
                    if (readingParam)
                    {
                        params = params + c;
                    }
                    else
                    {
                        blockName = blockName + c;
                    }
                }
            }

            if (blockName == slashStartTag)
            {
                break;
            }

            if ((blockName == compoundWord) && (params == kindFileStr))
            {
                // Keep entry
                word name;
                fileName path;
                word fName;
                bool foundName = false;
                bool foundPath = false;
                bool foundFName = false;
                bool earlyExit = false;
                while (!foundName || !foundPath || !foundFName)
                {
                    word entryName;
                    getEntry<word>(is, entryName);
                    if (entryName == nameWord)
                    {
                        getValue<word>(is, name);
                        if (nameRe.match(name))
                        {
                            foundName = true;
                        }
                        else
                        {
                            // Not interested in this compound
                            break;
                        }
                    }
                    else if (entryName == pathWord)
                    {
                        getValue<fileName>(is, path);

                        // Filter path on regExp
                        if (searchStrRe.match(path))
                        {
                            foundPath = true;
                        }
                        else
                        {
                            // Not interested in this compound
                            break;
                        }
                    }
                    else if (entryName == filenameWord)
                    {
                        getValue<word>(is, fName);
                        foundFName = true;
                    }
                    else
                    {
                        skipBlock(is, entryName);
                    }
                }

                if (foundPath && !earlyExit)
                {
                    word tName(path.components().last());

                    // Only insert if type is not already known
                    // NOTE: not ideal for cases where there are multiple types
                    //    but contained within different namespaces
                    //    preferentially take exact match if it exists
                    if (exactMatch && (tName + ".H") == name)
                    {
                        dictionary dict(dictionary::null);
                        dict.add("name", name);
                        dict.add("filename", fName + ".html");
                        dict.add("path", path);
                        this->add(tName, dict);
                    }
                    else if
                    (
                        !exactMatch
                     && !found(tName)
                     && wordRe(".*" + tName + ".*", wordRe::DETECT).match(name)
                    )
                    {
                        dictionary dict(dictionary::null);
                        dict.add("name", name);
                        dict.add("filename", fName + ".html");
                        dict.add("path", path);
                        this->add(tName, dict);
                    }
                }

                // Skip remanining entries
                skipBlock(is, blockName);
            }
            else
            {
                skipBlock(is, blockName);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::doxygenXmlParser::skipBlock
(
    IFstream& is,
    const word& blockName
) const
{
    // Recurse to move forward in 'is' until come across </blockName>
    string closeName;

    char c;
    while (is.good() && (closeName != blockName))
    {
        // Fast-forward until we reach a '<'
        while (is.get(c) && c  != '<')
        {}

        // Check to see if this is a closing block
        if (is.get(c) && c  == '/')
        {
            closeName = "";

            while (is.get(c) && c != '>')
            {
                closeName += c;
            }
        }
    }
}


void Foam::doxygenXmlParser::skipForward
(
    IFstream& is,
    const word& blockName
) const
{
    // Recurse to move forward in 'is' until come across <blockName>
    string entryName = "";
    char c;

    while (is.good() && (entryName != blockName))
    {
        entryName = "";

        // Fast-forward until we reach a '<'
        while (is.get(c) && c != '<')
        {}

        while (is.get(c) && c  != '>')
        {
            entryName = entryName + c;
        }
    }
}


// ************************************************************************* //
