/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "combustionModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::combustionModel> Foam::combustionModel::New
(
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
{
    IOobject combIO
    (
        IOobject
        (
            thermo.phasePropertyName(combustionProperties),
            thermo.T().mesh().time().constant(),
            thermo.T().mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word combModelName("none");
    if (combIO.typeHeaderOk<IOdictionary>(false))
    {
        IOdictionary(combIO).lookup("combustionModel") >> combModelName;
    }
    else
    {
        Info<< "Combustion model not active: "
            << thermo.phasePropertyName(combustionProperties)
            << " not found" << endl;
    }

    Info<< "Selecting combustion model " << combModelName << endl;

    const wordList cmpts2(basicThermo::splitThermoName(combModelName, 2));
    const wordList cmpts3(basicThermo::splitThermoName(combModelName, 3));
    if (cmpts2.size() == 2 || cmpts3.size() == 3)
    {
        combModelName = cmpts2.size() ? cmpts2[0] : cmpts3[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << combustionModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << combModelName << "." << endl;
    }

    const word thermoCombModelName =
        combModelName + '<' + thermo.thermoName() + '>';

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(combModelName);

    typename dictionaryConstructorTable::iterator thermoCstrIter =
        dictionaryConstructorTablePtr_->find(thermoCombModelName);

    if
    (
        cstrIter == dictionaryConstructorTablePtr_->end()
     && thermoCstrIter == dictionaryConstructorTablePtr_->end()
    )
    {
        FatalErrorInFunction
            << "Unknown " << combustionModel::typeName << " type "
            << combModelName << endl << endl;

        const wordList names(dictionaryConstructorTablePtr_->toc());

        wordList thisCmpts;
        thisCmpts.append(word::null);
        thisCmpts.append(basicThermo::splitThermoName(thermo.thermoName(), 5));

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
            << "Valid " << combustionModel::typeName << " types for this "
            << "thermodynamic model are:" << endl << validNames << endl;

        List<wordList> validCmpts1, validCmpts6;
        validCmpts1.append(wordList(1, word::null));
        validCmpts1[0][0] = combustionModel::typeName;
        validCmpts6.append(wordList(6, word::null));
        validCmpts6[0][0] = combustionModel::typeName;
        validCmpts6[0][1] = "transport";
        validCmpts6[0][2] = "thermo";
        validCmpts6[0][3] = "equationOfState";
        validCmpts6[0][4] = "specie";
        validCmpts6[0][5] = "energy";
        forAll(names, i)
        {
            const wordList cmpts1(basicThermo::splitThermoName(names[i], 1));
            const wordList cmpts6(basicThermo::splitThermoName(names[i], 6));
            if (cmpts1.size() == 1)
            {
                validCmpts1.append(cmpts1);
            }
            if (cmpts6.size() == 6)
            {
                validCmpts6.append(cmpts6);
            }
        }

        FatalErrorInFunction
            << "All " << validCmpts1[0][0]
            << " combinations are:" << endl << endl;
        printTable(validCmpts1, FatalErrorInFunction);

        FatalErrorInFunction << endl;

        FatalErrorInFunction
            << "All " << validCmpts6[0][0] << '/' << validCmpts6[0][1]
            << "/thermoPhysics combinations are:" << endl << endl;
        printTable(validCmpts6, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);
    }

    return autoPtr<combustionModel>
    (
        thermoCstrIter != dictionaryConstructorTablePtr_->end()
      ? thermoCstrIter()(combModelName, thermo, turb, combustionProperties)
      : cstrIter()(combModelName, thermo, turb, combustionProperties)
    );
}


// ************************************************************************* //
