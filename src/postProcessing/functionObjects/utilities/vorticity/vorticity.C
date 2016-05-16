/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "vorticity.H"
#include "volFields.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(vorticity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        vorticity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::vorticity::vorticity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    outputName_(typeName)
{
    read(dict);

    volVectorField* vorticityPtr
    (
        new volVectorField
        (
            IOobject
            (
                outputName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("0", dimless/dimTime, Zero)
        )
    );

    mesh_.objectRegistry::store(vorticityPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::vorticity::~vorticity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::vorticity::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");
    if (UName_ != "U")
    {
        outputName_ = typeName + "(" + UName_ + ")";
    }

    return true;
}


bool Foam::functionObjects::vorticity::execute(const bool postProcess)
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    volVectorField& vorticity = const_cast<volVectorField&>
    (
        mesh_.lookupObject<volVectorField>(outputName_)
    );

    vorticity = fvc::curl(U);

    return true;
}


bool Foam::functionObjects::vorticity::write(const bool postProcess)
{
    const volVectorField& vorticity =
        mesh_.lookupObject<volVectorField>(outputName_);

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << vorticity.name() << nl
        << endl;

    vorticity.write();

    return true;
}


// ************************************************************************* //
