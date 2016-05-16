/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "Q.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Q, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Q,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Q::Q
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    volScalarField* QPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless/sqr(dimTime), 0.0)
        )
    );

    mesh_.objectRegistry::store(QPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Q::~Q()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Q::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");

    return true;
}


bool Foam::functionObjects::Q::execute(const bool postProcess)
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);

    const volTensorField gradU(fvc::grad(U));

    volScalarField& Q =
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>(type())
        );

    Q = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));

    return true;
}


bool Foam::functionObjects::Q::write(const bool postProcess)
{
    const volScalarField& Q =
        mesh_.lookupObject<volScalarField>(type());

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << Q.name() << nl
        << endl;

    Q.write();

    return true;
}


// ************************************************************************* //
