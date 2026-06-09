/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "chemistryModel.H"
#include "basicThermo.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::chemistryModel> Foam::chemistryModel::New
(
    const fluidMulticomponentThermo& thermo
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.mesh().time().constant(),
            thermo.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word type
    (
        chemistryDict.lookupOrDefault<word>("type", "standard")
    );

    Info<< indentOrNl
        << "Selecting chemistryModel " << type << endl;

    // Check if the non-templated chemistryModel type exists on the table
    // providing support for third-party implementations which are independent
    // of OpenFOAM thermodynamics
    typename thermoConstructorTable::iterator cstrIter =
        thermoConstructorTablePtr_->find(type);

    if (cstrIter == thermoConstructorTablePtr_->end())
    {
        const word instantiatedType
        (
            type + '<' + thermo.thermoName() + ">"
        );

        Info<< "Instantiated as " << instantiatedType << endl;

        // Check if the chemistryModel type exists templated on
        // the current thermodynamics
        cstrIter = thermoConstructorTablePtr_->find(instantiatedType);

        if (cstrIter == thermoConstructorTablePtr_->end())
        {
            if
            (
                dynamicCode::allowSystemOperations
            && !dynamicCode::resolveTemplate
                (
                    chemistryModel::typeName
                ).empty()
            )
            {
                List<Pair<word>> substitutions
                (
                    basicThermo::thermoNameComponents(thermo.thermoName())
                );

                substitutions.append({"type", type});

                compileTemplate chemistryModel
                (
                    chemistryModel::typeName,
                    instantiatedType,
                    substitutions
                );
                cstrIter = thermoConstructorTablePtr_->find(instantiatedType);

                if (cstrIter == thermoConstructorTablePtr_->end())
                {
                    FatalIOErrorInFunction(chemistryDict)
                        << "Compilation and linkage of "
                        << chemistryModel::typeName << " type "
                        << type << nl << nl
                        << "failed." << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorInFunction(chemistryDict)
                    << "Unknown " << typeName_() << type
                    << " instantiated as " << instantiatedType
                    << endl;

                const wordList names(thermoConstructorTablePtr_->sortedToc());

                wordList thisCmpts(1, word::null);
                thisCmpts.append
                (
                    basicThermo::splitThermoName(thermo.thermoName(), 5)
                );

                List<wordList> validNames(1, {"type"});
                forAll(names, i)
                {
                    const wordList cmpts
                    (
                        basicThermo::splitThermoName(names[i], 6)
                    );

                    if
                    (
                        SubList<word>(cmpts, 5, 1)
                     == SubList<word>(thisCmpts, 5, 1)
                    )
                    {
                        validNames.append(SubList<word>(cmpts, 1));
                    }
                }

                FatalIOErrorInFunction(chemistryDict)
                    << "Valid chemistryModel types "
                       "for the current thermodynamics are:"
                    << nl << endl;
                printTable(validNames, FatalIOErrorInFunction(chemistryDict));

                FatalIOErrorInFunction(chemistryDict) << endl;

                List<wordList> validCmpts
                (
                    1,
                    {
                        "type",
                        "transport",
                        "thermo",
                        "equationOfState",
                        "specie",
                        "energy"
                    }
                );
                forAll(names, i)
                {
                    validCmpts.append
                    (
                        basicThermo::splitThermoName(names[i], 6)
                    );
                }

                FatalIOErrorInFunction(chemistryDict)
                    << "All chemistryModel/thermodynamics combinations are:"
                    << endl << endl;
                printTable(validCmpts, FatalIOErrorInFunction(chemistryDict));

                FatalIOErrorInFunction(chemistryDict) << exit(FatalIOError);
            }
        }
    }

    return autoPtr<chemistryModel>(cstrIter()(thermo));
}


// ************************************************************************* //
