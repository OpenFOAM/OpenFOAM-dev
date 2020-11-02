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

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    dictionary* dimensionedConstantsPtr_(nullptr);
}


Foam::dictionary& Foam::dimensionedConstants()
{
    return debug::switchSet
    (
        "DimensionedConstants",
        dimensionedConstantsPtr_
    );
}


void Foam::registerDimensionedConstant::lookup()
{
    dictionary& dict = dimensionedConstants();

    // Check that the entries exist.
    // Note: should make FatalError robust instead!

    if (!dict.found("unitSet"))
    {
        std::cerr<< "Cannot find unitSet in dictionary " << dict.name()
            << std::endl;
    }

    const word unitSetCoeffs(word(dict.lookup("unitSet")) + "Coeffs");

    if (!dict.found(unitSetCoeffs))
    {
        std::cerr<< "Cannot find " << unitSetCoeffs << " in dictionary "
            << dict.name() << std::endl;
    }

    dimensionedScalar::operator=
    (
        dimensionedScalar
        (
            name(),
            dimensions(),
            dict.subDict(unitSetCoeffs).subDict(group_)
        )
    );
}


void Foam::registerDimensionedConstantWithDefault::lookup()
{
    dictionary& dict = dimensionedConstants();

    const word unitSet(dict.lookup("unitSet"));
    dictionary& unitDict(dict.subDict(unitSet + "Coeffs"));

    const dimensionedScalar defaultValue(name(), defaultFunc_());

    if (unitDict.found(group_))
    {
        dictionary& groupDict = unitDict.subDict(group_);
        if (groupDict.found(name()))
        {
            dimensionedScalar::operator=
            (
                dimensionedScalar
                (
                    name(),
                    defaultValue.dimensions(),
                    groupDict.lookup(name())
                )
            );
        }
        else
        {
            groupDict.add(name(), defaultValue);
            dimensionedScalar::operator=(defaultValue);
        }
    }
    else
    {
        unitDict.add(group_, dictionary::null);
        unitDict.subDict(group_).add(name(), defaultValue);
        dimensionedScalar::operator=(defaultValue);
    }
}


void Foam::readDimensionedConstants(const dictionary& dict)
{
    if (dict.found("DimensionedConstants"))
    {
        InfoHeader
            << "Overriding DimensionedConstants according to "
            << dict.name() << endl;

        // Change dimensionedConstants dictionary in-memory
        dimensionedConstants().merge(dict.subDict("DimensionedConstants"));

        simpleObjectRegistry& objects = debug::dimensionedConstantObjects();

        IStringStream dummyIs("");

        forAllConstIter(simpleObjectRegistry, objects, iter)
        {
            const List<simpleRegIOobject*>& objects = *iter;

            forAll(objects, i)
            {
                objects[i]->readData(dummyIs);

                if (writeInfoHeader)
                {
                    Info<< "    ";
                    objects[i]->writeData(Info);
                    Info<< endl;
                }
            }
        }
    }
}


// ************************************************************************* //
