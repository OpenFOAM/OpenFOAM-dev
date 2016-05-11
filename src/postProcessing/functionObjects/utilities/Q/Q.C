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
#include "dictionary.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Q, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Q::Q
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

    volScalarField* QPtr
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

    mesh.objectRegistry::store(QPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Q::~Q()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::Q::read(const dictionary& dict)
{
    UName_ = dict.lookupOrDefault<word>("UName", "U");
}


void Foam::functionObjects::Q::execute()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const volVectorField& U =
        mesh.lookupObject<volVectorField>(UName_);

    const volTensorField gradU(fvc::grad(U));

    volScalarField& Q =
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(type())
        );

    Q = 0.5*(sqr(tr(gradU)) - tr(((gradU) & (gradU))));
}


void Foam::functionObjects::Q::end()
{
    execute();
}


void Foam::functionObjects::Q::timeSet()
{}


void Foam::functionObjects::Q::write()
{
    const volScalarField& Q =
        obr_.lookupObject<volScalarField>(type());

    Info<< type() << " " << name_ << " output:" << nl
        << "    writing field " << Q.name() << nl
        << endl;

    Q.write();
}


// ************************************************************************* //
