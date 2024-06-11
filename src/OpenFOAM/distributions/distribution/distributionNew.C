/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "distribution.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::distribution> Foam::distribution::New
(
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen,
    const bool report
)
{
    const word distributionType(dict.lookup("type"));

    if (report)
    {
        Info<< "Selecting " << typeName << " type " << distributionType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(distributionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type " << distributionType
            << nl << nl << "Valid " << typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<distribution> distributionPtr
    (
        cstrIter()
        (
            units,
            dict.optionalSubDict(distributionType + "Distribution"),
            sampleQ,
            std::move(rndGen)
        )
    );

    if (report)
    {
        Info<< incrIndent << indent
            << "min/average/max value = "
            << distributionPtr->min() << '/'
            << distributionPtr->mean() << '/'
            << distributionPtr->max()
            << decrIndent << endl;
    }

    return distributionPtr;
}


Foam::autoPtr<Foam::distribution> Foam::distribution::New
(
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    const randomGenerator::seed& s,
    const bool global,
    const bool report
)
{
    return New(units, dict, sampleQ, randomGenerator(s, global), report);
}


Foam::autoPtr<Foam::distribution> Foam::distribution::New
(
    autoPtr<distribution>& dPtr,
    const label sampleQ
)
{
    if (dPtr->sampleQ_ == sampleQ)
    {
        return autoPtr<distribution>(dPtr.ptr());
    }
    else
    {
        autoPtr<distribution> result(dPtr->clone(sampleQ));
        dPtr.clear();
        return result;
    }
}


// ************************************************************************* //
