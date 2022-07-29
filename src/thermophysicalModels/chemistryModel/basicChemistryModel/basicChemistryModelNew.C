/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
#include "compileTemplate.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicChemistryModel> Foam::basicChemistryModel::New
(
    const fluidMulticomponentThermo& thermo
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
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
            << "    chemistrySolver ode<chemistryModel<"
            << "rhoChemistryModel,sutherland<specie<janaf<perfectGas>,"
            << "sensibleInternalEnergy>>>>" << endl << endl << "becomes:"
            << endl << "    chemistryType" << endl << "    {" << endl
            << "        solver ode;" << endl << "    }" << exit(FatalError);
    }

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    const word solverName =
         chemistryTypeDict.lookupBackwardsCompatible<word>
         (
             {"solver", "chemistrySolver"}
         );

    const word methodName
    (
        chemistryTypeDict.lookupOrDefault<word>("method", "chemistryModel")
    );

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting chemistry solver " << chemistryTypeDictNew << endl;

    const word chemSolverNameName =
        solverName + '<' + methodName + '<' + thermo.thermoName() + ">>";

    typename thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(chemSolverNameName);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        if
        (
            dynamicCode::allowSystemOperations
         && !dynamicCode::resolveTemplate(basicChemistryModel::typeName).empty()
        )
        {
            List<Pair<word>> substitutions
            (
                basicThermo::thermoNameComponents(thermo.thermoName())
            );

            substitutions.append({"solver", solverName});
            substitutions.append({"method", methodName});

            compileTemplate chemistryModel
            (
                basicChemistryModel::typeName,
                chemSolverNameName,
                substitutions
            );
            cstrIter = thermoConstructorTablePtr_->find(chemSolverNameName);

            if (cstrIter == thermoConstructorTablePtr_->end())
            {
                FatalErrorInFunction
                    << "Compilation and linkage of "
                    << basicChemistryModel::typeName << " type " << nl
                    << "chemistryType" << chemistryTypeDict << nl << nl
                    << "failed." << exit(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown " << typeName_() << " type " << chemistryTypeDictNew
                << endl;

            const wordList names(thermoConstructorTablePtr_->sortedToc());

            wordList thisCmpts;
            thisCmpts.append(word::null);
            thisCmpts.append(word::null);
            thisCmpts.append
            (
                basicThermo::splitThermoName(thermo.thermoName(), 5)
            );

            List<wordList> validNames;
            validNames.append(wordList(2, word::null));
            validNames[0][0] = "solver";
            validNames[0][1] = "method";
            forAll(names, i)
            {
                const wordList cmpts(basicThermo::splitThermoName(names[i], 7));

                if
                (
                    SubList<word>(cmpts, 5, 2)
                 == SubList<word>(thisCmpts, 5, 2)
                )
                {
                    validNames.append(SubList<word>(cmpts, 2));
                }
            }

            FatalErrorInFunction
                << "Valid " << validNames[0][0] << '/' << validNames[0][1]
                << " combinations for this thermodynamic model are:"
                << endl << endl;
            printTable(validNames, FatalErrorInFunction);

            FatalErrorInFunction << endl;

            List<wordList> validCmpts;
            validCmpts.append(wordList(7, word::null));
            validCmpts[0][0] = "solver";
            validCmpts[0][1] = "method";
            validCmpts[0][2] = "transport";
            validCmpts[0][3] = "thermo";
            validCmpts[0][4] = "equationOfState";
            validCmpts[0][5] = "specie";
            validCmpts[0][6] = "energy";
            forAll(names, i)
            {
                validCmpts.append(basicThermo::splitThermoName(names[i], 7));
            }

            FatalErrorInFunction
                << "All " << validCmpts[0][0] << '/' << validCmpts[0][1]
                << "/thermodynamics combinations are:"
                << endl << endl;
            printTable(validCmpts, FatalErrorInFunction);

            FatalErrorInFunction << exit(FatalError);
        }
    }

    return autoPtr<basicChemistryModel>(cstrIter()(thermo));
}


// ************************************************************************* //
