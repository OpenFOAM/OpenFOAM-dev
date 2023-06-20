/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    foamToC

Description
    Run-time selection table of contents printing and interrogation.

    The run-time selection tables are populated by the optionally specified
    solver class and any additional libraries listed in the \c -libs option or
    all libraries using the \c -allLibs option.  Once populated the tables can
    be searched and printed by a range of options listed below.  Table entries
    are printed with the corresponding library they are in to aid selection
    and the addition of \c libs entries to ensure availability to the solver.

Usage
    \b foamToC [OPTION]
      - \par -noLibs
        Load only the core libOpenFOAM.so library by default

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Pre-load additional libraries

      - \par -solvers \<name\>
        List solvers

      - \par -solver \<name\>
        Load libraries associated with specified solver

      - \par -listLibs
        List libraries as they are loaded

      - \par switches,
        List all available debug, info and optimisation switches

      - \par all,
        List the contents of all the run-time selection tables

      - \par tables
        List the run-time selection table names (this is the default action)

      - \par table \<name\>
        List the contents of the specified table or the list sub-tables

      - \par search \<name\> or \"\<regular expression\>\"
        Search for and list the tables containing the given name
        or all matches to the given regular expression

      - \par scalarBCs,
        List scalar field boundary conditions (fvPatchField<scalar>)

      - \par vectorBCs,
        List vector field boundary conditions (fvPatchField<vector>)

      - \par functionObjects,
        List functionObjects

      - \par fvModels,
        List fvModels

      - \par fvConstraints,
        List fvConstraints

    Example usage:
      - Print the list of solvers
        \verbatim
            foamToC -solvers
        \endverbatim

      - Print the list of scalar boundary conditions (fvPatchField<scalar>)
        provided by the \c fluid solver without additional libraries:
        \verbatim
            foamToC -solver fluid -scalarBCs
        \endverbatim

      - Print the list of RAS momentum transport models provided by the
        \c fluid solver:
        \verbatim
            foamToC -solver fluid -table RAScompressibleMomentumTransportModel
        \endverbatim

      - Print the list of functionObjects provided by the
        \c multicomponentFluid solver with the libfieldFunctionObjects.so
        library:
        \verbatim
            foamToC -solver multicomponentFluid \
                -libs '("libfieldFunctionObjects.so")' -functionObjects
        \endverbatim

      - Print a complete list of all run-time selection tables:
        \verbatim
            foamToC -tables
            or simply
            foamToC
        \endverbatim

      - Print a complete list of all entries in all run-time selection tables:
        \verbatim
            foamToC -all
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "runTimeSelectionToC.H"
#include "dlLibraryTable.H"
#include "HashSet.H"
#include "IOmanip.H"

using namespace Foam;

// Enable run-time selection table of contents caching
// Overrides the enableRunTimeSelectionToC = false in libOpenFOAM
bool Foam::debug::enableRunTimeSelectionToC = true;


HashTable<HashTable<wordHashSet>> baseTypeNameToC()
{
    HashTable<HashTable<wordHashSet>> baseTypeNameToC;

    forAllConstIter
    (
        debug::runTimeSelectionToCType,
        debug::runTimeSelectionToC,
        iter
    )
    {
        const word& baseType = iter.key();
        const word& baseTypeName = iter().first();

        if (!baseTypeNameToC.found(baseTypeName))
        {
            baseTypeNameToC.insert
            (
                baseTypeName,
                HashTable<wordHashSet>()
            );
        }

        if (!baseTypeNameToC[baseTypeName].found(baseType))
        {
            baseTypeNameToC[baseTypeName].insert
            (
                baseType,
                iter().second()
            );
        }
    }

    return baseTypeNameToC;
}


void printToC(const word& tableName)
{
    bool found = false;

    if (debug::runTimeSelectionToC.found(tableName))
    {
        const wordList toc
        (
            debug::runTimeSelectionToC[tableName].second().sortedToc()
        );

        Info<< "Contents of table " << tableName;

        if (debug::runTimeSelectionToC[tableName].first() != tableName)
        {
            Info<< ", base type "
                << debug::runTimeSelectionToC[tableName].first();
        }

        Info<< ":" << endl;

        forAll(toc, i)
        {
            Info<< "    " << setf(ios_base::left) << setw(40) << toc[i]
                << debug::runTimeSelectionToC[tableName].second()[toc[i]]
                << endl;
        }

        found = true;
    }
    else
    {
        const HashTable<HashTable<wordHashSet>> runTimeSelectionToC
        (
            baseTypeNameToC()
        );

        if (runTimeSelectionToC.found(tableName))
        {
            const wordList toc(runTimeSelectionToC[tableName].sortedToc());

            Info<< "Tables of type " << tableName << endl;
            forAll(toc, i)
            {
                Info<< "    " << toc[i] << endl;
            }

            found = true;
        }
    }

    if (!found)
    {
        Info<< "Table " << tableName << " not found" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    entry::disableFunctionEntries = true;
    writeInfoHeader = false;

    argList::initValidTables::clear();

    argList::addBoolOption
    (
        "noLibs",
        "Load only the core libOpenFOAM.so library by default"
    );

    argList::addOption
    (
        "libs",
        "'(\"lib1.so\" ... \"libN.so\")'",
        "Pre-load additional libraries"
    );

    argList::addBoolOption
    (
        "solvers",
        "List solvers"
    );

    argList::addOption
    (
        "solver",
        "name",
        "Load libraries associated with specified solver"
    );

    argList::addBoolOption
    (
        "listLibs",
        "List libraries as they are loaded"
    );

    argList::addBoolOption
    (
        "switches",
        "List all available debug, info and optimisation switches"
    );

    argList::addBoolOption
    (
        "all",
        "List the contents of all the run-time selection tables"
    );

    argList::addBoolOption
    (
        "tables",
        "List the run-time selection table names"
    );

    argList::addOption
    (
        "table",
        "name",
        "List the contents of the specified table"
    );

    argList::addOption
    (
        "search",
        "name",
        "List the tables containing the given name"
    );

    argList::addBoolOption
    (
        "scalarBCs",
        "List scalar field boundary conditions (fvPatchField<scalar>)"
    );

    argList::addBoolOption
    (
        "vectorBCs",
        "List vector field boundary conditions (fvPatchField<vector>)"
    );

    argList::addBoolOption
    (
        "functionObjects",
        "List functionObjects"
    );

    argList::addBoolOption
    (
        "fvModels",
        "List fvModels"
    );

    argList::addBoolOption
    (
        "fvConstraints",
        "List fvConstraints"
    );

    const argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    word solverName;
    if (args.optionReadIfPresent("solver", solverName))
    {
        libs.open("lib" + solverName + ".so");
    }

    const string libDir(getEnv("FOAM_LIBBIN"));
    const fileNameList libNames(readDir(libDir));

    const bool listLibs = args.optionFound("listLibs");

    if (!args.optionFound("noLibs") && solverName == word::null)
    {
        if (listLibs)
        {
            Info << "Loading libraries:" << nl;
        }
        forAll(libNames, i)
        {
            if (libNames[i].ext() == "so")
            {
                if (listLibs)
                {
                    Info << "    " << libNames[i].c_str() << nl;
                }
                libs.open(libDir/libNames[i]);
            }
        }
        Info << endl;
    }

    bool done = false;

    if (args.optionFound("switches"))
    {
        debug::listSwitches();
        done = true;
    }

    word tableName;
    if (args.optionReadIfPresent("table", tableName))
    {
        printToC(tableName);
        done = true;
    }

    word name;
    if (args.optionReadIfPresent("search", name))
    {
        wordList names(name);

        // If the name is a regular expression pattern
        // search for all matching occurrences
        if (wordRe::isPattern(name))
        {
            wordRe nameRe(name);
            nameRe.compile();

            HashSet<word> unsortedNames;

            forAllConstIter
            (
                debug::runTimeSelectionToCType,
                debug::runTimeSelectionToC,
                iter
            )
            {
                forAllConstIter(HashTable<word>, iter().second(), iter2)
                {
                    if (nameRe.match(iter2.key()))
                    {
                        unsortedNames.insert(iter2.key());
                    }
                }
            }

            names = unsortedNames.toc();
        }

        forAll(names, i)
        {
            const word& name = names[i];

            HashTable<HashTable<word>> baseTypeNameTables;

            forAllConstIter
            (
                debug::runTimeSelectionToCType,
                debug::runTimeSelectionToC,
                iter
            )
            {
                const word& baseType = iter.key();
                const word& baseTypeName = iter().first();

                if (iter().second().found(name))
                {
                    if (!baseTypeNameTables.found(baseTypeName))
                    {
                        baseTypeNameTables.insert
                        (
                            baseTypeName,
                            HashTable<word>()
                        );
                    }

                    baseTypeNameTables[baseTypeName].insert
                    (
                        baseType,
                        iter().second()[name]
                    );
                }
            }

            if (baseTypeNameTables.size())
            {
                const wordList toc(baseTypeNameTables.sortedToc());

                if (toc.size() == 1)
                {
                    Info<< name << " is in table " << endl;
                }
                else if (toc.size() > 1)
                {
                    Info<< name << " is in tables " << endl;
                }

                forAll(toc, i)
                {
                    if
                    (
                        baseTypeNameTables[toc[i]].size() == 1
                     && toc[i] == baseTypeNameTables[toc[i]].begin().key()
                    )
                    {
                        Info<< "    " << setf(ios_base::left) << setw(40)
                            << toc[i] << baseTypeNameTables[toc[i]].begin()()
                            << endl;
                    }
                    else
                    {
                        const wordList tocc
                        (
                            baseTypeNameTables[toc[i]].sortedToc()
                        );

                        Info<< "    " << toc[i] << endl;
                        forAll(tocc, j)
                        {
                            Info<< "        "
                                << setf(ios_base::left) << setw(40) << tocc[j]
                                << baseTypeNameTables[toc[i]][tocc[j]]
                                << endl;
                        }
                    }
                }
            }
            else
            {
                Info<< name << " not found" << endl;
            }
        }

        done = true;
    }

    if (args.optionFound("all"))
    {
        Info<< "ToC:" << nl
            << debug::runTimeSelectionToC << endl;
        done = true;
    }

    if (args.optionFound("solvers"))
    {
        printToC("solver");
        done = true;
    }

    if (args.optionFound("scalarBCs"))
    {
        Info<< "Scalar boundary conditions:" << endl;
        printToC("fvPatchScalarField");
        done = true;
    }

    if (args.optionFound("vectorBCs"))
    {
        Info<< "vector boundary conditions:" << endl;
        printToC("fvPatchVectorField");
        done = true;
    }

    if (args.optionFound("functionObjects"))
    {
        Info<< "functionObjects:" << endl;
        printToC("functionObject");
        done = true;
    }

    if (args.optionFound("fvModels"))
    {
        Info<< "fvModels:" << endl;
        printToC("fvModel");
        done = true;
    }

    if (args.optionFound("fvConstraints"))
    {
        Info<< "fvConstraints:" << endl;
        printToC("fvConstraint");
        done = true;
    }

    if (args.optionFound("tables") || !done)
    {
        const HashTable<HashTable<wordHashSet>> runTimeSelectionToC
        (
            baseTypeNameToC()
        );

        const wordList toc(runTimeSelectionToC.sortedToc());

        Info<< "Run-time selection tables:" << nl;

        forAll(toc, i)
        {
            if
            (
                runTimeSelectionToC[toc[i]].size() == 1
             && toc[i] == runTimeSelectionToC[toc[i]].begin().key()
            )
            {
                Info<< "    " << toc[i] << endl;
            }
            else
            {
                const wordList tocc
                (
                    runTimeSelectionToC[toc[i]].sortedToc()
                );

                Info<< "    " << toc[i] << endl;
                forAll(tocc, j)
                {
                    Info<< "        " << tocc[j] << endl;
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
