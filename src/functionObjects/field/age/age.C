/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "age.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "inletOutletFvPatchField.H"
#include "wallFvPatch.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(age, 0);
    addToRunTimeSelectionTable(functionObject, age, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::age::patchTypes() const
{
    wordList result
    (
        mesh_.boundary().size(),
        inletOutletFvPatchField<scalar>::typeName
    );

    forAll(mesh_.boundary(), patchi)
    {
        if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
        {
            result[patchi] = zeroGradientFvPatchField<scalar>::typeName;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::age::age
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phiName_(),
    rhoName_(),
    nCorr_(0),
    schemesField_()

{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::age::~age()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::age::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    dict.readIfPresent("nCorr", nCorr_);

    schemesField_ = dict.lookupOrDefault<word>("schemesField", typeName);

    return true;
}


bool Foam::functionObjects::age::execute()
{
    return true;
}


bool Foam::functionObjects::age::write()
{
    volScalarField t
    (
        IOobject
        (
            typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, 0),
        patchTypes()
    );

    const word divScheme("div(phi," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    // This only works because the null constructed inletValue for an
    // inletOutletFvPatchField is zero. If we needed any other value we would
    // have to loop over the inletOutlet patches and explicitly set the
    // inletValues. We would need to change the interface of inletOutlet in
    // order to do this.

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (label i = 0; i <= nCorr_; ++ i)
        {
            fvScalarMatrix tEqn
            (
                fvm::div(phi, t, divScheme) == rho
            );

            tEqn.relax(relaxCoeff);

            tEqn.solve(schemesField_);
        }
    }
    else
    {
        for (label i = 0; i <= nCorr_; ++ i)
        {
            fvScalarMatrix tEqn
            (
                fvm::div(phi, t, divScheme) == dimensionedScalar(1)
            );

            tEqn.relax(relaxCoeff);

            tEqn.solve(schemesField_);
        }
    }

    Info<< "Min/max age:" << min(t).value() << ' ' << max(t).value() << endl;

    t.write();

    return true;
}


// ************************************************************************* //
