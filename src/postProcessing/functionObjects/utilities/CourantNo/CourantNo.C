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

#include "CourantNo.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CourantNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        CourantNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::functionObjects::CourantNo::byRho
(
    const tmp<volScalarField::Internal>& Co
) const
{
    if (Co().dimensions() == dimDensity)
    {
        return Co/obr_.lookupObject<volScalarField>(rhoName_);
    }
    else
    {
        return Co;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::CourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    phiName_("phi"),
    rhoName_("rho")
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    volScalarField* CourantNoPtr
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
            dimensionedScalar("0", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    mesh.objectRegistry::store(CourantNoPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::~CourantNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CourantNo::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");

    return true;
}


bool Foam::functionObjects::CourantNo::execute(const bool postProcess)
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    volScalarField& Co = const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>(type())
    );

    Co.ref() = byRho
    (
        (0.5*mesh.time().deltaT())
       *fvc::surfaceSum(mag(phi))()()
       /mesh.V()
    );
    Co.correctBoundaryConditions();

    return true;
}


bool Foam::functionObjects::CourantNo::write(const bool postProcess)
{
    const volScalarField& CourantNo =
        obr_.lookupObject<volScalarField>(type());

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << CourantNo.name() << nl
        << endl;

    CourantNo.write();

    return true;
}


// ************************************************************************* //
