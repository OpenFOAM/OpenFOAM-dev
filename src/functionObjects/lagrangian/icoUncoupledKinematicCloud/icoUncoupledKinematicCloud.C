/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "icoUncoupledKinematicCloud.H"
#include "singlePhaseTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(icoUncoupledKinematicCloud, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        icoUncoupledKinematicCloud,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::icoUncoupledKinematicCloud::icoUncoupledKinematicCloud
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    g_
    (
        IOobject
        (
            "g",
            time_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    laminarTransport_
    (
        mesh_.lookupObject<singlePhaseTransportModel>("transportProperties")
    ),
    rhoValue_
    (
        "rho",
        dimDensity,
        laminarTransport_
    ),
    rho_
    (
        IOobject
        (
            "rho",
            time_.timeName(),
            mesh_
        ),
        mesh_,
        rhoValue_
    ),
    mu_("mu", rhoValue_*laminarTransport_.nu()),
    U_
    (
        mesh_.lookupObject<volVectorField>(dict.lookupOrDefault<word>("U", "U"))
    ),
    kinematicCloudName_
    (
        dict.lookupOrDefault<word>("kinematicCloud", "kinematicCloud")
    ),
    kinematicCloud_
    (
        kinematicCloudName_,
        rho_,
        U_,
        mu_,
        g_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::icoUncoupledKinematicCloud::~icoUncoupledKinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::icoUncoupledKinematicCloud::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::icoUncoupledKinematicCloud::execute()
{
    mu_ = rhoValue_*laminarTransport_.nu();
    kinematicCloud_.evolve();

    return true;
}


bool Foam::functionObjects::icoUncoupledKinematicCloud::write()
{
    return true;
}


// ************************************************************************* //
