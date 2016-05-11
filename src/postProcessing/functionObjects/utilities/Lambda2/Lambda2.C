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
#include "dictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Lambda2, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Lambda2::Lambda2
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    UName_("U")
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField* Lambda2Ptr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimless/sqr(dimTime), 0.0)
        )
    );

    mesh.objectRegistry::store(Lambda2Ptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Lambda2::~Lambda2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::Lambda2::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");
}


void Foam::functionObjects::Lambda2::execute()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const volVectorField& U =
        mesh.lookupObject<volVectorField>(UName_);

    const volTensorField gradU(fvc::grad(U));

    const volTensorField SSplusWW
    (
        (symm(gradU) & symm(gradU))
      + (skew(gradU) & skew(gradU))
    );

    volScalarField& Lambda2 =
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(type())
        );

    Lambda2 = -eigenValues(SSplusWW)().component(vector::Y);
}


void Foam::functionObjects::Lambda2::end()
{
    execute();
}


void Foam::functionObjects::Lambda2::timeSet()
{}


void Foam::functionObjects::Lambda2::write()
{
    const volScalarField& Lambda2 =
        obr_.lookupObject<volScalarField>(type());

    Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << Lambda2.name() << nl
        << endl;

    Lambda2.write();
}


// ************************************************************************* //
