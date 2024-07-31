/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "constantTemperature.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::saturationTemperatureModel>
Foam::saturationTemperatureModel::New(const dictionary& dict)
{
    return New(NullObjectRef<word>(), dict);
}


Foam::autoPtr<Foam::saturationTemperatureModel>
Foam::saturationTemperatureModel::New
(
    const word& name,
    const dictionary& dict
)
{
    if (!isNull(name) && !dict.isDict(name))
    {
        Istream& is(dict.lookup(name, false));

        token firstToken(is);
        if (!firstToken.isWord())
        {
            return autoPtr<saturationTemperatureModel>
            (
                new saturationModels::constantTemperature
                (
                    dimensionedScalar(name, dimTemperature, dict)
                )
            );
        }
    }

    const bool isType = isNull(name);
    const bool isDict = !isType && dict.isDict(name);

    const word modelTypeName =
        isType ? dict.lookup("type")
      : isDict ? dict.subDict(name).lookup("type")
      : dict.lookup<word>(name);

    const dictionary& coeffDict =
        isType ? dict
      : isDict ? dict.subDict(name)
      : dict.optionalSubDict(name + "Coeffs");

    Info<< "Selecting " << typeName << " " << modelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " << type "
            << modelTypeName << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(coeffDict);
}


// ************************************************************************* //
