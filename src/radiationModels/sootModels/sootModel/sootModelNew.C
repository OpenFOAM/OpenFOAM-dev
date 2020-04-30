/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "error.H"
#include "sootModel.H"
#include "noSoot.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModels::sootModel>
Foam::radiationModels::sootModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    // Get the soot model type name
    word modelType(sootModels::noSoot::typeName);
    if (dict.found(sootModel::typeName))
    {
        dict.lookup(sootModel::typeName) >> modelType;
        Info<< "Selecting soot model " << modelType << endl;
    }
    const wordList cmpts(basicThermo::splitThermoName(modelType, 3));
    if (cmpts.size() == 3)
    {
        modelType = cmpts[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << sootModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << modelType << "." << endl;
    }

    // Get the thermo model type names
    word thermoType(word::null);
    if (mesh.foundObject<basicThermo>(basicThermo::dictName))
    {
        const basicThermo& thermo =
            mesh.lookupObject<basicThermo>(basicThermo::dictName);

        thermoType = thermo.thermoName();
    }

    // Construct a thermo-soot model type name
    const word thermoModelType = modelType + '<' + thermoType + '>';

    // Lookup both possible model names
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    dictionaryConstructorTable::iterator thermoCstrIter =
        dictionaryConstructorTablePtr_->find(thermoModelType);

    // Construct and return
    if (thermoCstrIter != dictionaryConstructorTablePtr_->end())
    {
        return autoPtr<sootModel>(thermoCstrIter()(dict, mesh, modelType));
    }
    else if (cstrIter != dictionaryConstructorTablePtr_->end())
    {
        return autoPtr<sootModel>(cstrIter()(dict, mesh, modelType));
    }
    else
    {
        FatalErrorInFunction
            << "Unknown " << sootModel::typeName << " type "
            << modelType << nl << nl;

        const wordList names(dictionaryConstructorTablePtr_->sortedToc());

        wordList thisCmpts;
        thisCmpts.append(word::null);
        thisCmpts.append(basicThermo::splitThermoName(thermoType, 5));

        wordList validNames;
        forAll(names, i)
        {
            wordList cmpts(basicThermo::splitThermoName(names[i], 1));
            if (cmpts.size() != 1)
            {
                cmpts = basicThermo::splitThermoName(names[i], 6);
            }

            bool isValid = true;
            for (label i = 1; i < cmpts.size() && isValid; ++ i)
            {
                isValid = isValid && cmpts[i] == thisCmpts[i];
            }

            if (isValid)
            {
                validNames.append(cmpts[0]);
            }
        }

        FatalErrorInFunction
            << "Valid " << sootModel::typeName << " types for this "
            << "thermodynamic model are:" << endl << validNames << endl;

        List<wordList> validCmpts;
        validCmpts.append(wordList(6, word::null));
        validCmpts[0][0] = sootModel::typeName;
        validCmpts[0][1] = "transport";
        validCmpts[0][2] = "thermo";
        validCmpts[0][3] = "equationOfState";
        validCmpts[0][4] = "specie";
        validCmpts[0][5] = "energy";
        forAll(names, i)
        {
            const wordList cmpts1(basicThermo::splitThermoName(names[i], 1));
            const wordList cmpts6(basicThermo::splitThermoName(names[i], 6));
            if (cmpts1.size() == 1)
            {
                validCmpts.append(wordList(6, "<any>"));
                validCmpts.last()[0] = cmpts1[0];
            }
            if (cmpts6.size() == 6)
            {
                validCmpts.append(cmpts6);
            }
        }

        FatalErrorInFunction
            << "All " << sootModel::typeName
            << "/thermoPhysics combinations are:" << endl << endl;
        printTable(validCmpts, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);

       return autoPtr<sootModel>(nullptr);
    }
}


// ************************************************************************* //
