/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "XiReactionRate.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(XiReactionRate, 0);
    addToRunTimeSelectionTable(functionObject, XiReactionRate, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::XiReactionRate::XiReactionRate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::XiReactionRate::~XiReactionRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::XiReactionRate::fields() const
{
    return wordList{"b", "Su", "Xi"};
}


bool Foam::functionObjects::XiReactionRate::execute()
{
    return true;
}


bool Foam::functionObjects::XiReactionRate::write()
{
    const volScalarField& b =
        mesh_.lookupObject<volScalarField>("b");

    const volScalarField& Su =
        mesh_.lookupObject<volScalarField>("Su");

    const volScalarField& Xi =
        mesh_.lookupObject<volScalarField>("Xi");

    volScalarField St
    (
        IOobject
        (
            "St",
            time_.name(),
            mesh_
        ),
        Xi*Su
    );

    Log << "    Writing turbulent flame-speed field " << St.name()
        << " to " << time_.name() << endl;

    St.write();

    volScalarField wdot
    (
        IOobject
        (
            "wdot",
            time_.name(),
            mesh_
        ),
        St*mag(fvc::grad(b))
    );

    Log << "    Writing reaction-rate field " << wdot.name()
        << " to " << time_.name() << endl;

    wdot.write();

    return true;
}


// ************************************************************************* //
