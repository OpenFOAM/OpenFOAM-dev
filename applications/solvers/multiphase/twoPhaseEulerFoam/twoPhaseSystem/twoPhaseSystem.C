/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "BlendedInterfacialModel.H"
#include "virtualMassModel.H"
#include "heatTransferModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "fvMatrix.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fixedValueFvsPatchFields.H"
#include "blendingMethod.H"
#include "HashPtrTable.H"
#include "UniformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phase1_
    (
        *this,
        *this,
        wordList(lookup("phases"))[0]
    ),

    phase2_
    (
        *this,
        *this,
        wordList(lookup("phases"))[1]
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->calcPhi()
    ),

    dgdt_
    (
        IOobject
        (
            "dgdt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dgdt", dimless/dimTime, 0)
    )
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);


    // Blending
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(lookup("phases"))
            )
        );
    }


    // Pairs

    phasePair::scalarTable sigmaTable(lookup("sigma"));
    phasePair::dictTable aspectRatioTable(lookup("aspectRatio"));

    pair_.set
    (
        new phasePair
        (
            phase1_,
            phase2_,
            g,
            sigmaTable
        )
    );

    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_,
            phase2_,
            g,
            sigmaTable,
            aspectRatioTable
        )
    );

    pair2In1_.set
    (
        new orderedPhasePair
        (
            phase2_,
            phase1_,
            g,
            sigmaTable,
            aspectRatioTable
        )
    );


    // Models

    drag_.set
    (
        new BlendedInterfacialModel<dragModel>
        (
            lookup("drag"),
            (
                blendingMethods_.found("drag")
              ? blendingMethods_["drag"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_,
            false // Do not zero drag coefficient at fixed-flux BCs
        )
    );

    virtualMass_.set
    (
        new BlendedInterfacialModel<virtualMassModel>
        (
            lookup("virtualMass"),
            (
                blendingMethods_.found("virtualMass")
              ? blendingMethods_["virtualMass"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    heatTransfer_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            lookup("heatTransfer"),
            (
                blendingMethods_.found("heatTransfer")
              ? blendingMethods_["heatTransfer"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    lift_.set
    (
        new BlendedInterfacialModel<liftModel>
        (
            lookup("lift"),
            (
                blendingMethods_.found("lift")
              ? blendingMethods_["lift"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    wallLubrication_.set
    (
        new BlendedInterfacialModel<wallLubricationModel>
        (
            lookup("wallLubrication"),
            (
                blendingMethods_.found("wallLubrication")
              ? blendingMethods_["wallLubrication"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    turbulentDispersion_.set
    (
        new BlendedInterfacialModel<turbulentDispersionModel>
        (
            lookup("turbulentDispersion"),
            (
                blendingMethods_.found("turbulentDispersion")
              ? blendingMethods_["turbulentDispersion"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::rho() const
{
    return phase1_*phase1_.thermo().rho() + phase2_*phase2_.thermo().rho();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::U() const
{
    return phase1_*phase1_.U() + phase2_*phase2_.U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_)*phase1_.phi()
      + fvc::interpolate(phase2_)*phase2_.phi();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kd() const
{
    return drag_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Kdf() const
{
    return drag_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Vm() const
{
    return virtualMass_->K();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Vmf() const
{
    return virtualMass_->Kf();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::Kh() const
{
    return heatTransfer_->K();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseSystem::F() const
{
    return lift_->F<vector>() + wallLubrication_->F<vector>();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseSystem::Ff() const
{
    return lift_->Ff() + wallLubrication_->Ff();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseSystem::D() const
{
    return turbulentDispersion_->D();
}


void Foam::twoPhaseSystem::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();


    surfaceScalarField phic("phic", phi_);
    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alpha1alpha2f;

    if (pPrimeByA_.valid())
    {
        alpha1alpha2f =
            fvc::interpolate(max(alpha1, scalar(0)))
           *fvc::interpolate(max(alpha2, scalar(0)));

        surfaceScalarField phiP
        (
            pPrimeByA_()*fvc::snGrad(alpha1, "bounded")*mesh_.magSf()
        );

        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dgdt_.dimensions(), 0.0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        forAll(dgdt_, celli)
        {
            if (dgdt_[celli] > 0.0)
            {
                Sp[celli] -= dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
                Su[celli] += dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
            }
            else if (dgdt_[celli] < 0.0)
            {
                Sp[celli] += dgdt_[celli]/max(alpha1[celli], 1e-4);
            }
        }

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        phase1_.correctInflowOutflow(alphaPhic1);

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    UniformField<scalar>(phase1_.alphaMax()),
                    zeroField()
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_.alphaPhi() = alphaPhic10;
                }
                else
                {
                    phase1_.alphaPhi() += alphaPhic10;
                }
            }

            phase1_.alphaPhi() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                Sp,
                Su,
                UniformField<scalar>(phase1_.alphaMax()),
                zeroField()
            );

            phase1_.alphaPhi() = alphaPhic1;
        }

        if (pPrimeByA_.valid())
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1alpha2f()*pPrimeByA_(), alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_.alphaPhi() += alpha1Eqn.flux();
        }

        phase1_.alphaRhoPhi() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhi() = phi_ - phase1_.alphaPhi();
        phase2_.correctInflowOutflow(phase2_.alphaPhi());
        phase2_.alphaRhoPhi() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.max(0);
        alpha1.min(1);

        alpha2 = scalar(1) - alpha1;
    }
}


void Foam::twoPhaseSystem::correct()
{
    phase1_.correct();
    phase2_.correct();
}


void Foam::twoPhaseSystem::correctTurbulence()
{
    phase1_.turbulence().correct();
    phase2_.turbulence().correct();
}


bool Foam::twoPhaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_.read(*this);
        readOK &= phase2_.read(*this);

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


const Foam::dimensionedScalar& Foam::twoPhaseSystem::sigma() const
{
    return pair_->sigma();
}


// ************************************************************************* //
