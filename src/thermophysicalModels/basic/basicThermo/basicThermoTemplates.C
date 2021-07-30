/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "basicThermo.H"
#include "wordIOList.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupCstrIter
(
    const dictionary& thermoTypeDict,
    Table* tablePtr,
    const int nCmpt,
    const char* cmptNames[],
    const word& thermoTypeName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);

    if (cstrIter == tablePtr->end())
    {
        if
        (
            dynamicCode::allowSystemOperations
         && !dynamicCode::resolveTemplate(Thermo::typeName).empty()
        )
        {
            compileTemplate thermo
            (
                Thermo::typeName,
                thermoTypeName,
                List<Pair<word>>
                {
                    {"type", thermoTypeDict.lookup("type")},
                    {"mixture", thermoTypeDict.lookup("mixture")},
                    {"transport", thermoTypeDict.lookup("transport")},
                    {"thermo", thermoTypeDict.lookup("thermo")},
                    {
                        "equationOfState",
                        thermoTypeDict.lookup("equationOfState")
                    },
                    {"specie", thermoTypeDict.lookup("specie")},
                    {"energy", thermoTypeDict.lookup("energy")}
                }
            );
            cstrIter = tablePtr->find(thermoTypeName);

            if (cstrIter == tablePtr->end())
            {
                FatalErrorInFunction
                    << "Compilation and linkage of "
                    << Thermo::typeName << " type " << nl
                    << "thermoType" << thermoTypeDict << nl << nl
                    << "failed." << nl << nl
                    << "Valid " << Thermo::typeName << " types are:"
                    << nl << nl;
            }
        }
        else
        {
            // Print error message if package not found in the table
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " type " << nl
                << "thermoType" << thermoTypeDict << nl << nl
                << "Valid " << Thermo::typeName << " types are:"
                << nl << nl;
        }

        if (cstrIter == tablePtr->end())
        {
            // Get the list of all the suitable thermo packages available
            wordList validThermoTypeNames(tablePtr->sortedToc());

            // Build a table of the thermo packages constituent parts
            DynamicList<wordList> validThermoTypeNameCmpts;

            // Set row zero to the column headers
            validThermoTypeNameCmpts.append(wordList(nCmpt));
            forAll(validThermoTypeNameCmpts[0], i)
            {
                validThermoTypeNameCmpts[0][i] = cmptNames[i];
            }

            // Split the thermo package names into their constituent parts and
            // add them to the table, removing any incompatible entries from the
            // list
            forAll(validThermoTypeNames, i)
            {
                const wordList names
                (
                    Thermo::splitThermoName(validThermoTypeNames[i], nCmpt)
                );

                if (names.size())
                {
                    validThermoTypeNameCmpts.append(names);
                }
            }

            // Print the table of available packages
            printTable(validThermoTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }
    }

    return cstrIter;
}


template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupCstrIter
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        if (thermoTypeDict.found("properties"))
        {
            const int nCmpt = 4;
            const char* cmptNames[nCmpt] =
            {
                "type",
                "mixture",
                "properties",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                word(thermoTypeDict.lookup("type")) + '<'
              + word(thermoTypeDict.lookup("mixture")) + '<'
              + word(thermoTypeDict.lookup("properties")) + ','
              + word(thermoTypeDict.lookup("energy")) + ">>"
            );

            return lookupCstrIter<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                nCmpt,
                cmptNames,
                thermoTypeName
            );
        }
        else
        {
            const int nCmpt = 7;
            const char* cmptNames[nCmpt] =
            {
                "type",
                "mixture",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                word(thermoTypeDict.lookup("type")) + '<'
              + word(thermoTypeDict.lookup("mixture")) + '<'
              + word(thermoTypeDict.lookup("transport")) + '<'
              + word(thermoTypeDict.lookup("thermo")) + '<'
              + word(thermoTypeDict.lookup("equationOfState")) + '<'
              + word(thermoTypeDict.lookup("specie")) + ">>,"
              + word(thermoTypeDict.lookup("energy")) + ">>>"
            );

            return lookupCstrIter<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                nCmpt,
                cmptNames,
                thermoTypeName
            );
        }
    }
    else
    {
        const word thermoTypeName(thermoDict.lookup("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);

        if (cstrIter == tablePtr->end())
        {
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " type "
                << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl
                << tablePtr->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    const IOdictionary thermoDict
    (
        physicalProperties::findModelDict(mesh, phaseName)
    );

    typename Thermo::fvMeshConstructorTable::iterator cstrIter =
        lookupCstrIter<Thermo, typename Thermo::fvMeshConstructorTable>
        (
            thermoDict,
            Thermo::fvMeshConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, phaseName));
}


// ************************************************************************* //
