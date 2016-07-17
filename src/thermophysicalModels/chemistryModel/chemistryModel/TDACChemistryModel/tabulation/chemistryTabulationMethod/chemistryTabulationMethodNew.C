/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "chemistryTabulationMethod.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistryTabulationMethod<CompType, ThermoType>>
Foam::chemistryTabulationMethod<CompType, ThermoType>::New
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
{
    IOdictionary thermoDict
    (
         IOobject
         (
              "thermophysicalProperties",
              dict.db().time().constant(),
              dict.db(),
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
        FatalIOErrorInFunction(thermoDict)
            << "thermoType is in the old format and must be upgraded"
            << exit(FatalIOError);
    }

    dictionary tabdict(dict.subDict("tabulation"));

    word chemistryTabulationMethodName =
        word(tabdict.lookup("method")) + '<'
      + word(dict.subDict("chemistryType").lookup("chemistryThermo")) + ','
      + thermoTypeName + '>';

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistryTabulationMethodName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown chemistryTabulationMethodType type "
            << chemistryTabulationMethodName
            << endl << endl
            << "Valid chemistryTabulationMethodType types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<chemistryTabulationMethod<CompType, ThermoType>>
    (
        cstrIter()(dict, chemistry)
    );
}


// ************************************************************************* //
