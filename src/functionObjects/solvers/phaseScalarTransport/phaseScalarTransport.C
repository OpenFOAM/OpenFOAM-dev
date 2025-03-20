/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "fixedValueFvPatchField.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "nonOrthogonalSolutionControl.H"
#include "phaseScalarTransport.H"
#include "surfaceFields.H"
#include "momentumTransportModel.H"
#include "wallFvPatch.H"
#include "zeroGradientFvPatchField.H"

#define PhiDimensionErrorInFunction(phi)                                       \
    FatalErrorInFunction                                                       \
        << "Incompatible dimensions for " << phi.name() << ": "                \
        << phi.dimensions() << nl                                              \
        << "Dimensions should be " << dimMass/dimTime << " or "                \
        << dimVolume/dimTime << exit(FatalError)

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phaseScalarTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        phaseScalarTransport,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField& Foam::functionObjects::phaseScalarTransport::Phi()
{
    if (!PhiPtr_.valid())
    {
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);
        const volScalarField& p =
            mesh_.lookupObject<volScalarField>(pName_);

        wordList PhiPatchFieldTypes(mesh_.boundaryMesh().size());
        forAll(p.boundaryField(), patchi)
        {
            PhiPatchFieldTypes[patchi] =
                p.boundaryField()[patchi].fixesValue()
              ? fixedValueFvPatchField<scalar>::typeName
              : zeroGradientFvPatchField<scalar>::typeName;
        }

        PhiPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    typedName(IOobject::groupName("Phi", phaseName_)),
                    time_.name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(phi.dimensions()/dimLength, Zero),
                PhiPatchFieldTypes
            )
        );

        mesh_.schemes().setFluxRequired(PhiPtr_->name());
    }

    return PhiPtr_();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::phaseScalarTransport::alphaPhi()
{
    if (!solveAlphaPhi_)
    {
        return mesh_.lookupObject<surfaceScalarField>(alphaPhiName_);
    }

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Make a crude guess of the phase flux using default interpolation
    tmp<surfaceScalarField> tAlphaPhi
    (
        new surfaceScalarField
        (
            alphaPhiName_,
            phi*fvc::interpolate(alpha)
        )
    );
    surfaceScalarField& alphaPhi = tAlphaPhi.ref();

    // Get the potential field
    volScalarField& Phi(this->Phi());

    // Construct the scheme names
    const word laplacianScheme = "laplacian(" + pName_ + ")";

    // Debug writing. Write the material derivative of alpha, before and after
    // the solution of the potential and the correction of alphaPhi. Before
    // correction the field should be non-zero, and after it should be
    // comparable to the solution tolerance.
    auto writeDDt = [&](const label i)
    {
        const volScalarField DDtAlpha
        (
            "DDt("
          + IOobject::groupName
            (
                IOobject::member(alpha.name()) + Foam::name(i),
                IOobject::group(alpha.name())
            )
          + ")",
            fvc::ddt(alpha) + fvc::div(alphaPhi)
        );
        Info<< type() << ": Writing " << DDtAlpha.name() << endl;
        DDtAlpha.write();
    };
    if (debug && time_.writeTime())
    {
        writeDDt(0);
    }

    // Lookup the non-orthogonal solution control
    nonOrthogonalSolutionControl& control =
        mesh_.lookupObjectRef<nonOrthogonalSolutionControl>
        (
            solutionControl::typeName
        );

    // Solve for the potential and correct alphaPhi with the resulting flux
    if (phi.dimensions() == dimVolume/dimTime)
    {
        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(Phi, laplacianScheme)
              + fvc::ddt(alpha)
              + fvc::div(alphaPhi)
            );

            PhiEqn.solve(pName_);

            if (control.finalNonOrthogonalIter())
            {
                alphaPhi += PhiEqn.flux();
            }
        }
    }
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        while (control.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(Phi, laplacianScheme)
              + fvc::ddt(rho, alpha)
              + fvc::div(alphaPhi)
            );

            PhiEqn.solve(pName_);

            if (control.finalNonOrthogonalIter())
            {
                alphaPhi += PhiEqn.flux();
            }
        }
    }
    else
    {
        PhiDimensionErrorInFunction(phi);
    }

    // Debug writing
    if (debug && time_.writeTime())
    {
        writeDDt(1);
    }

    return tAlphaPhi;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::phaseScalarTransport::D
(
    const surfaceScalarField& alphaPhi
) const
{
    const word Dname("D" + s_.name());

    if (diffusivity_ == scalarTransport::diffusivityType::constant)
    {
        return volScalarField::New
        (
            Dname,
            mesh_,
            dimensionedScalar(dimKinematicViscosity, D_)
        );
    }
    else
    {
        const word& nameNoPhase = momentumTransportModel::typeName;
        const word namePhase = IOobject::groupName(nameNoPhase, phaseName_);

        // Try looking up the phase transport model, then try the mixture
        // transport model, then fail with an error relating to the phase
        // transport model
        const momentumTransportModel& turbulence =
            mesh_.foundObject<momentumTransportModel>(namePhase)
          ? mesh_.lookupObject<momentumTransportModel>(namePhase)
          : mesh_.foundObject<momentumTransportModel>(nameNoPhase)
          ? mesh_.lookupObject<momentumTransportModel>(nameNoPhase)
          : mesh_.lookupObject<momentumTransportModel>(namePhase);

        return volScalarField::New
        (
            Dname,
            alphal_*turbulence.nu() + alphat_*turbulence.nut()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phaseScalarTransport::phaseScalarTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookup("field")),
    phaseName_(IOobject::group(fieldName_)),
    s_
    (
        IOobject
        (
            fieldName_,
            time_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    PhiPtr_(nullptr)
{
    if (phaseName_ == word::null)
    {
        FatalErrorInFunction
            << "Field \"" << fieldName_ << "\" does not have a phase extension "
            << "in its name. If it is associated with \"phaseA\" then it "
            << "should be named \"" << fieldName_ << ".phaseA\"."
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::phaseScalarTransport::~phaseScalarTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phaseScalarTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    solveAlphaPhi_ = dict.lookupOrDefault<bool>("solveAlphaPhi", false);

    alphaName_ =
        dict.lookupOrDefault<word>
        (
            "alpha",
            IOobject::groupName("alpha", phaseName_)
        );
    const word defaultAlphaPhiName =
        IOobject::groupName("alphaPhi", phaseName_);
    alphaPhiName_ =
        solveAlphaPhi_
      ? typedName(defaultAlphaPhiName)
      : dict.lookupOrDefault<word>("alphaPhi", defaultAlphaPhiName);
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ =
        dict.lookupOrDefault<word>
        (
            "rho",
            IOobject::groupName("rho", phaseName_)
        );
    pName_ = dict.lookupOrDefault<word>("p", "p");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_);

    diffusivity_ =
        scalarTransport::diffusivityTypeNames_.read(dict.lookup("diffusivity"));

    switch(diffusivity_)
    {
        case scalarTransport::diffusivityType::none:
            break;

        case scalarTransport::diffusivityType::constant:
            dict.lookup("D") >> D_;
            break;

        case scalarTransport::diffusivityType::viscosity:
            alphal_ =
                dict.lookupBackwardsCompatible<scalar>({"alphal", "alphaD"});
            alphat_ =
                dict.lookupBackwardsCompatible<scalar>({"alphat", "alphaDt"});
            break;
    }

    nCorr_ = dict.lookupOrDefault<label>("nCorr", 0);
    residualAlpha_ = dict.lookupOrDefault<scalar>("residualAlpha", rootSmall);
    writeAlphaField_ = dict.lookupOrDefault<bool>("writeAlphaField", true);

    return true;
}


Foam::wordList Foam::functionObjects::phaseScalarTransport::fields() const
{
    return wordList{alphaName_, alphaPhiName_, phiName_, pName_};
}


bool Foam::functionObjects::phaseScalarTransport::execute()
{
    Info<< type() << ": Executing" << endl;

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    // Get the phase flux
    tmp<surfaceScalarField> tAlphaPhi(this->alphaPhi());
    const surfaceScalarField& alphaPhi = tAlphaPhi();

    // Get the relaxation coefficient
    const scalar relaxCoeff =
        mesh_.solution().relaxEquation(schemesField_)
      ? mesh_.solution().equationRelaxationFactor(schemesField_)
      : 0;

    // Models and constraints
    const Foam::fvModels& fvModels = Foam::fvModels::New(mesh_);
    const Foam::fvConstraints& fvConstraints = Foam::fvConstraints::New(mesh_);

    // Solve
    if (alphaPhi.dimensions() == dimVolume/dimTime)
    {
        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix fieldEqn
            (
                fvm::ddt(alpha, s_)
              + fvm::div
                (
                    alphaPhi,
                    s_,
                    "div(" + alphaPhi.name() + "," + schemesField_ + ")"
                )
             ==
                fvModels.source(alpha, s_)
              - fvm::ddt(residualAlpha_, s_)
              + fvc::ddt(residualAlpha_, s_)
            );

            if (diffusivity_ != scalarTransport::diffusivityType::none)
            {
                const volScalarField D(this->D(alphaPhi));

                fieldEqn -=
                    fvm::laplacian
                    (
                        fvc::interpolate(alpha)*fvc::interpolate(D),
                        s_,
                        "laplacian(" + D.name() + "," + schemesField_ + ")"
                    );
            }

            fieldEqn.relax(relaxCoeff);
            fvConstraints.constrain(fieldEqn);
            fieldEqn.solve(schemesField_);
            fvConstraints.constrain(s_);
        }
    }
    else if (alphaPhi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix fieldEqn
            (
                fvm::ddt(alpha, rho, s_)
              + fvm::div
                (
                    alphaPhi,
                    s_,
                    "div(" + alphaPhi.name() + "," + schemesField_ + ")"
                )
             ==
                fvModels.source(alpha, rho, s_)
              - fvm::ddt(residualAlpha_*rho, s_)
              + fvc::ddt(residualAlpha_*rho, s_)
            );

            if (diffusivity_ != scalarTransport::diffusivityType::none)
            {
                const volScalarField D(this->D(alphaPhi));

                fieldEqn -=
                    fvm::laplacian
                    (
                        fvc::interpolate(alpha)*fvc::interpolate(rho*D),
                        s_,
                        "laplacian(" + D.name() + "," + schemesField_ + ")"
                    );
            }

            fieldEqn.relax(relaxCoeff);
            fvConstraints.constrain(fieldEqn);
            fieldEqn.solve(schemesField_);
            fvConstraints.constrain(s_);
        }
    }
    else
    {
        PhiDimensionErrorInFunction(alphaPhi);
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::phaseScalarTransport::write()
{
    s_.write();

    if (writeAlphaField_)
    {
        const volScalarField& alpha =
            mesh_.lookupObject<volScalarField>(alphaName_);

        volScalarField alphaS
        (
            IOobject
            (
                "alpha" + fieldName_.capitalise(),
                time_.name(),
                mesh_
            ),
            alpha*s_
        );

        alphaS.write();
    }

    if (PhiPtr_.valid())
    {
        PhiPtr_->write();
    }

    return true;
}


// ************************************************************************* //
