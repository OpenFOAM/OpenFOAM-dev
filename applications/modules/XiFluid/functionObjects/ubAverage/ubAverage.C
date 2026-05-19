/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "ubAverage.H"
#include "XiFluid.H"
#include "ubRhoThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ubAverage, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ubAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ubAverage::calc()
{
    if (foundObject<solvers::XiFluid>(solver::typeName))
    {
        const solvers::XiFluid& XiFluid
        (
            lookupObject<solvers::XiFluid>(solver::typeName)
        );

        const ubRhoThermo& ubThermo = XiFluid.thermo;

        const word fuName
        (
            IOobject::groupName(fieldName_, ubRhoThermo::unburntPhaseName)
        );

        const word fbName
        (
            IOobject::groupName(fieldName_, ubRhoThermo::burntPhaseName)
        );

        if
        (
            foundObject<volScalarField>(fuName)
         && foundObject<volScalarField>(fbName)
        )
        {
            const volScalarField& fu = lookupObject<volScalarField>(fuName);
            const volScalarField& fb = lookupObject<volScalarField>(fbName);

            store
            (
                resultName_,
                ubThermo.alphau()*fu + ubThermo.alphab()*fb
            );

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ubAverage::ubAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ubAverage::~ubAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::ubAverage::fields() const
{
    return
    {
        IOobject::groupName(fieldName_, ubRhoThermo::unburntPhaseName),
        IOobject::groupName(fieldName_, ubRhoThermo::burntPhaseName)
    };
}


// ************************************************************************* //
