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

#include "Lambda2.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Lambda2, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Lambda2,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Lambda2::Lambda2
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);

    volScalarField* Lambda2Ptr
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

    mesh_.objectRegistry::store(Lambda2Ptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Lambda2::~Lambda2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::Lambda2::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");

    return true;
}


bool Foam::functionObjects::Lambda2::execute(const bool postProcess)
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);

    const volTensorField gradU(fvc::grad(U));

    const volTensorField SSplusWW
    (
        (symm(gradU) & symm(gradU))
      + (skew(gradU) & skew(gradU))
    );

    volScalarField& Lambda2 =
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>(type())
        );

    Lambda2 = -eigenValues(SSplusWW)().component(vector::Y);

    return true;
}


bool Foam::functionObjects::Lambda2::write(const bool postProcess)
{
    const volScalarField& Lambda2 =
        obr_.lookupObject<volScalarField>(type());

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << Lambda2.name() << nl
        << endl;

    Lambda2.write();

    return true;
}


// ************************************************************************* //
