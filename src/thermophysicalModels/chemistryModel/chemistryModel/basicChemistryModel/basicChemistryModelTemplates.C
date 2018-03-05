/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    typename ChemistryModel::reactionThermo& thermo
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!chemistryDict.isDict("chemistryType"))
    {
        FatalErrorInFunction
            << "Template parameter based chemistry solver selection is no "
            << "longer supported. Please create a chemistryType dictionary"
            << "instead." << endl << endl << "For example, the entry:" << endl
            << "    chemistrySolver ode<StandardChemistryModel<"
            << "rhoChemistryModel,sutherlandspecie<janaf<perfectGas>,"
            << "sensibleInternalEnergy>>" << endl << endl << "becomes:" << endl
            << "    chemistryType" << endl << "    {" << endl
            << "        solver ode;" << endl << "        method standard;"
            << endl << "    }" << exit(FatalError);
    }

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    const word& solverName
    (
        chemistryTypeDict.found("solver")
      ? chemistryTypeDict.lookup("solver")
      : chemistryTypeDict.found("chemistrySolver")
      ? chemistryTypeDict.lookup("chemistrySolver")
      : chemistryTypeDict.lookup("solver") // error if neither entry is found
    );

    const word& methodName
    (
        chemistryTypeDict.lookupOrDefault<word>
        (
            "method",
            chemistryTypeDict.lookupOrDefault<bool>("TDAC", false)
          ? "TDAC"
          : "standard"
        )
    );

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting chemistry solver " << chemistryTypeDictNew << endl;

    typedef typename ChemistryModel::thermoConstructorTable cstrTableType;
    cstrTableType* cstrTable = ChemistryModel::thermoConstructorTablePtr_;

    const word chemSolverCompThermoName =
        solverName + '<' + methodName + '<'
      + ChemistryModel::reactionThermo::typeName + ','
      + thermo.thermoName() + ">>";

    typename cstrTableType::iterator cstrIter =
        cstrTable->find(chemSolverCompThermoName);

    if (cstrIter == cstrTable->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName_() << " type " << solverName << endl
            << endl;

        const wordList names(cstrTable->toc());

        wordList thisCmpts;
        thisCmpts.append(word::null);
        thisCmpts.append(word::null);
        thisCmpts.append(ChemistryModel::reactionThermo::typeName);
        thisCmpts.append(basicThermo::splitThermoName(thermo.thermoName(), 5));

        List<wordList> validNames;
        validNames.append(wordList(2, word::null));
        validNames[0][0] = "solver";
        validNames[0][1] = "method";
        forAll(names, i)
        {
            const wordList cmpts(basicThermo::splitThermoName(names[i], 8));

            bool isValid = true;
            for (label i = 2; i < cmpts.size() && isValid; ++ i)
            {
                isValid = isValid && cmpts[i] == thisCmpts[i];
            }

            if (isValid)
            {
                validNames.append(SubList<word>(cmpts, 2));
            }
        }

        FatalErrorInFunction
            << "All " << validNames[0][0] << '/' << validNames[0][1]
            << "combinations for this thermodynamic model are:"
            << endl << endl;
        printTable(validNames, FatalErrorInFunction);

        FatalErrorInFunction << endl;

        List<wordList> validCmpts;
        validCmpts.append(wordList(8, word::null));
        validCmpts[0][0] = "solver";
        validCmpts[0][1] = "method";
        validCmpts[0][2] = "reactionThermo";
        validCmpts[0][3] = "transport";
        validCmpts[0][4] = "thermo";
        validCmpts[0][5] = "equationOfState";
        validCmpts[0][6] = "specie";
        validCmpts[0][7] = "energy";
        forAll(names, i)
        {
            validCmpts.append(basicThermo::splitThermoName(names[i], 8));
        }

        FatalErrorInFunction
            << "All " << validCmpts[0][0] << '/' << validCmpts[0][1] << '/'
            << validCmpts[0][2] << "/thermoPhysics combinations are:"
            << endl << endl;
        printTable(validCmpts, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);
    }

    return autoPtr<ChemistryModel>(cstrIter()(thermo));
}

// ************************************************************************* //
