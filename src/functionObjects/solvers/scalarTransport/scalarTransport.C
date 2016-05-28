/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "scalarTransport.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fvScalarMatrix.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(scalarTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        scalarTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::scalarTransport::boundaryTypes() const
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    wordList bTypes(U.boundaryField().size());

    forAll(bTypes, patchi)
    {
        const fvPatchField<vector>& pf = U.boundaryField()[patchi];
        if (isA<fixedValueFvPatchVectorField>(pf))
        {
            bTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
        else
        {
            bTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
        }
    }

    return bTypes;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::scalarTransport::DT
(
    const surfaceScalarField& phi
) const
{
    typedef incompressible::turbulenceModel icoModel;
    typedef compressible::turbulenceModel cmpModel;

    if (userDT_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "DT",
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DT", phi.dimensions()/dimLength, DT_)
            )
        );
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
    {
        const icoModel& model = mesh_.lookupObject<icoModel>
        (
            turbulenceModel::propertiesName
        );

        return model.nuEff();
    }
    else if (mesh_.foundObject<cmpModel>(turbulenceModel::propertiesName))
    {
        const cmpModel& model = mesh_.lookupObject<cmpModel>
        (
            turbulenceModel::propertiesName
        );

        return model.muEff();
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "DT",
                    mesh_.time().timeName(),
                    mesh_.time(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("DT", phi.dimensions()/dimLength, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransport::scalarTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    DT_(0.0),
    nCorr_(0),
    fvOptions_(mesh_),
    T_
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0),
        boundaryTypes()
    )
{
    read(dict);

    if (resetOnStartUp_)
    {
        T_ == dimensionedScalar("zero", dimless, 0.0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::scalarTransport::~scalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::scalarTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    UName_ = dict.lookupOrDefault<word>("U", "U");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    userDT_ = false;
    if (dict.readIfPresent("DT", DT_))
    {
        userDT_ = true;
    }

    dict.lookup("resetOnStartUp") >> resetOnStartUp_;

    dict.readIfPresent("nCorr", nCorr_);

    dict.lookup("autoSchemes") >> autoSchemes_;

    fvOptions_.reset(dict.subDict("fvOptions"));

    return true;
}


bool Foam::functionObjects::scalarTransport::execute(const bool postProcess)
{
    Info<< type() << " output:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // calculate the diffusivity
    volScalarField DT(this->DT(phi));

    // set schemes
    word schemeVar = T_.name();
    if (autoSchemes_)
    {
        schemeVar = UName_;
    }

    word divScheme("div(phi," + schemeVar + ")");
    word laplacianScheme("laplacian(" + DT.name() + "," + schemeVar + ")");

    // set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemeVar))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemeVar);
    }

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        // solve
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(rho, T_)
              + fvm::div(phi, T_, divScheme)
              - fvm::laplacian(DT, T_, laplacianScheme)
             ==
                fvOptions_(rho, T_)
            );

            TEqn.relax(relaxCoeff);

            fvOptions_.constrain(TEqn);

            TEqn.solve(mesh_.solverDict(schemeVar));
        }
    }
    else if (phi.dimensions() == dimVolume/dimTime)
    {
        // solve
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T_)
              + fvm::div(phi, T_, divScheme)
              - fvm::laplacian(DT, T_, laplacianScheme)
             ==
                fvOptions_(T_)
            );

            TEqn.relax(relaxCoeff);

            fvOptions_.constrain(TEqn);

            TEqn.solve(mesh_.solverDict(schemeVar));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << endl;
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::scalarTransport::write(const bool postProcess)
{
    return true;
}


// ************************************************************************* //
