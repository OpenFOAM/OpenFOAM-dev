/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "momentumTransportModel.H"
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


bool Foam::functionObjects::age::converged
(
    const int nCorr,
    const scalar initialResidual
) const
{
    if (initialResidual < tolerance_)
    {
        Info<< "Field " << typeName
            << " converged in " << nCorr << " correctors\n" << endl;

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::age::age
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

Foam::functionObjects::age::~age()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::age::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    nCorr_ = dict.lookupOrDefault<int>("nCorr", 5);
    schemesField_ = dict.lookupOrDefault<word>("schemesField", typeName);
    diffusion_ = dict.lookupOrDefault<Switch>("diffusion", false);
    tolerance_ = dict.lookupOrDefault<scalar>("tolerance", 1e-5);

    return true;
}


Foam::wordList Foam::functionObjects::age::fields() const
{
    return wordList{phiName_, rhoName_};
}


bool Foam::functionObjects::age::execute()
{
    tmp<volScalarField> tage
    (
        new volScalarField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(dimTime, 0),
            patchTypes()
        )
    );
    volScalarField& age = tage.ref();

    const word divScheme("div(phi," + schemesField_ + ")");

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.solution().relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.solution().equationRelaxationFactor(schemesField_);
    }

    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(mesh_)
    );

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

        const word laplacianScheme("laplacian(muEff," + schemesField_ + ")");

        tmp<volScalarField> tnuEff;
        if (diffusion_)
        {
            tnuEff =
                mesh_.lookupObject<momentumTransportModel>
                (
                    momentumTransportModel::typeName
                ).nuEff();
        }

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix ageEqn
            (
                fvm::div(phi, age, divScheme) == rho + fvModels.source(rho, age)
            );

            if (diffusion_)
            {
                ageEqn -= fvm::laplacian(rho*tnuEff(), age, laplacianScheme);
            }

            ageEqn.relax(relaxCoeff);

            fvConstraints.constrain(ageEqn);

            if (converged(i, ageEqn.solve(schemesField_).initialResidual()))
            {
                break;
            };

            fvConstraints.constrain(age);
        }
    }
    else
    {
        tmp<volScalarField> tnuEff;
        word laplacianScheme;

        if (diffusion_)
        {
            tnuEff =
                mesh_.lookupObject<momentumTransportModel>
                (
                    momentumTransportModel::typeName
                ).nuEff();

            laplacianScheme =
                "laplacian(" + tnuEff().name() + ',' + schemesField_ + ")";
        }

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix ageEqn
            (
                fvm::div(phi, age, divScheme)
             == dimensionedScalar(1) + fvModels.source(age)
            );

            if (diffusion_)
            {
                ageEqn -= fvm::laplacian(tnuEff(), age, laplacianScheme);
            }

            ageEqn.relax(relaxCoeff);

            fvConstraints.constrain(ageEqn);

            if (converged(i, ageEqn.solve(schemesField_).initialResidual()))
            {
                break;
            }

            fvConstraints.constrain(age);
        }
    }

    Info<< "Min/max age:" << min(age).value()
        << ' ' << max(age).value() << endl;

    return store(tage);
}


bool Foam::functionObjects::age::write()
{
    return writeObject(typeName);
}


// ************************************************************************* //
