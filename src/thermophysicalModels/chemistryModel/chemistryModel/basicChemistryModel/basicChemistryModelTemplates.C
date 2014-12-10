/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "basicChemistryModel.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::autoPtr<ChemistryModel> Foam::basicChemistryModel::New
(
    const fvMesh& mesh
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word chemistryTypeName;

    if (chemistryDict.isDict("chemistryType"))
    {
        const dictionary& chemistryTypeDict
        (
            chemistryDict.subDict("chemistryType")
        );

        Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

        const int nCmpt = 7;
        const char* cmptNames[nCmpt] =
        {
            "chemistrySolver",
            "chemistryThermo",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        };

        IOdictionary thermoDict
        (
            IOobject
            (
                "thermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        word thermoTypeName;

        if (thermoDict.isDict("thermoType"))
        {
            const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
            thermoTypeName =
                word(thermoTypeDict.lookup("transport")) + '<'
              + word(thermoTypeDict.lookup("thermo")) + '<'
              + word(thermoTypeDict.lookup("equationOfState")) + '<'
              + word(thermoTypeDict.lookup("specie")) + ">>,"
              + word(thermoTypeDict.lookup("energy")) + ">";
        }
        else
        {
             FatalIOErrorIn
             (
                 (ChemistryModel::typeName + "::New(const mesh&)").c_str(),
                 thermoDict
             )   << "thermoType is in the old format and must be upgraded"
                 << exit(FatalIOError);
        }

        // Construct the name of the chemistry type from the components
        chemistryTypeName =
            word(chemistryTypeDict.lookup("chemistrySolver")) + '<'
          + word(chemistryTypeDict.lookup("chemistryThermo")) + ','
          + thermoTypeName + ">";

        typename ChemistryModel::fvMeshConstructorTable::iterator cstrIter =
            ChemistryModel::fvMeshConstructorTablePtr_->find(chemistryTypeName);

        if (cstrIter == ChemistryModel::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorIn(ChemistryModel::typeName + "::New(const mesh&)")
                << "Unknown " << ChemistryModel::typeName << " type " << nl
                << "chemistryType" << chemistryTypeDict << nl << nl
                << "Valid " << ChemistryModel ::typeName << " types are:"
                << nl << nl;

            // Get the list of all the suitable chemistry packages available
            wordList validChemistryTypeNames
            (
                ChemistryModel::fvMeshConstructorTablePtr_->sortedToc()
            );

            // Build a table of the thermo packages constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validChemistryTypeNameCmpts
            (
                validChemistryTypeNames.size() + 1
            );

            validChemistryTypeNameCmpts[0].setSize(nCmpt);
            forAll(validChemistryTypeNameCmpts[0], j)
            {
                validChemistryTypeNameCmpts[0][j] = cmptNames[j];
            }

            // Split the thermo package names into their constituent parts
            forAll(validChemistryTypeNames, i)
            {
                validChemistryTypeNameCmpts[i+1] = basicThermo::splitThermoName
                (
                    validChemistryTypeNames[i],
                    nCmpt
                );
            }

            // Print the table of available packages
            // in terms of their constituent parts
            printTable(validChemistryTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }

        return autoPtr<ChemistryModel>(cstrIter()(mesh));
    }
    else
    {
        chemistryTypeName =
            word(chemistryDict.lookup("chemistryType"));

        Info<< "Selecting chemistry type " << chemistryTypeName << endl;

        typename ChemistryModel::fvMeshConstructorTable::iterator cstrIter =
            ChemistryModel::fvMeshConstructorTablePtr_->find(chemistryTypeName);

        if (cstrIter == ChemistryModel::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorIn(ChemistryModel::typeName + "::New(const mesh&)")
                << "Unknown " << ChemistryModel::typeName << " type "
                << chemistryTypeName << nl << nl
                << "Valid ChemistryModel types are:" << nl
                << ChemistryModel::fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return autoPtr<ChemistryModel>(cstrIter()(mesh));
    }
}

// ************************************************************************* //
