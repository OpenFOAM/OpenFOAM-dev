/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "momentumTransferSystem.H"

#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "generateBlendedInterfacialModels.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "fvcMeshPhi.H"
#include "fvcReconstruct.H"

#include "scalarMatrices.H"
#include "pimpleNoLoopControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentumTransferSystem, 0);
}


const Foam::word Foam::momentumTransferSystem::propertiesName
(
    "momentumTransfer"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::IOobject Foam::momentumTransferSystem::io(const phaseSystem& fluid) const
{
    typeIOobject<IOdictionary> result
    (
        propertiesName,
        fluid.mesh().time().constant(),
        fluid.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    // Check for and warn about momentum transfer models being in the old
    // location in constant/phaseProperties
    const wordList modelEntries
    ({
        modelName<blendedDragModel>(),
        modelName<blendedVirtualMassModel>(),
        modelName<blendedLiftModel>(),
        modelName<blendedWallLubricationModel>(),
        modelName<blendedTurbulentDispersionModel>()
    });

    OStringStream modelEntriesString;
    forAll(modelEntries, i)
    {
        modelEntriesString<< modelEntries[i];
        if (i < modelEntries.size() - 2) modelEntriesString<< ", ";
        if (i == modelEntries.size() - 2) modelEntriesString<< " and ";
    }

    if (!result.headerOk())
    {
        result.readOpt() = IOobject::NO_READ;

        WarningInFunction
            << "Specifying momentum transfer model entries - "
            << modelEntriesString.str().c_str() << " - in "
            << fluid.relativeObjectPath() << " is deprecated. These "
            << "entries should now be specified in "
            << result.relativeObjectPath() << "." << endl;
    }
    else
    {
        forAll(modelEntries, i)
        {
            if (!fluid.found(modelEntries[i])) continue;

            WarningInFunction
                << "Momentum transfer model entries - "
                << modelEntriesString.str().c_str() << " - in "
                << fluid.relativeObjectPath() << " are no longer used. These "
                << "entries are now read from " << result.relativeObjectPath()
                << "." << endl;

            break;
        }
    }

    return result;
}


template<class ModelType>
const Foam::dictionary& Foam::momentumTransferSystem::modelsDict() const
{
    const word key = modelName<ModelType>();

    return
        readOpt() == IOobject::NO_READ && fluid_.found(key)
      ? fluid_.subDict(key)
      : subDict(key);
}


void Foam::momentumTransferSystem::addTmpField
(
    tmp<surfaceScalarField>& result,
    const tmp<surfaceScalarField>& field
) const
{
    if (result.valid())
    {
        result.ref() += field;
    }
    else
    {
        if (field.isTmp())
        {
            result = field;
        }
        else
        {
            result = field().clone();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentumTransferSystem::momentumTransferSystem(const phaseSystem& fluid)
:
    IOdictionary(io(fluid)),
    fluid_(fluid),
    dragModels_
    (
        generateBlendedInterfacialModels<blendedDragModel>
        (
            fluid_,
            modelsDict<blendedDragModel>()
        )
    ),
    virtualMassModels_
    (
        generateBlendedInterfacialModels<blendedVirtualMassModel>
        (
            fluid_,
            modelsDict<blendedVirtualMassModel>()
        )
    ),
    liftModels_
    (
        generateBlendedInterfacialModels<blendedLiftModel>
        (
            fluid_,
            modelsDict<blendedLiftModel>()
        )
    ),
    wallLubricationModels_
    (
        generateBlendedInterfacialModels<blendedWallLubricationModel>
        (
            fluid_,
            modelsDict<blendedWallLubricationModel>()
        )
    ),
    turbulentDispersionModels_
    (
        generateBlendedInterfacialModels<blendedTurbulentDispersionModel>
        (
            fluid_,
            modelsDict<blendedTurbulentDispersionModel>()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentumTransferSystem::~momentumTransferSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::PtrList<Foam::surfaceScalarField>
Foam::momentumTransferSystem::Fs() const
{
    PtrList<surfaceScalarField> Fs(fluid_.phases().size());

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const phaseInterface& interface = liftModelIter()->interface();

        const volVectorField F(liftModelIter()->F());

        addField
        (
            interface.phase1(),
            "F",
            fvc::flux(F),
            Fs
        );
        addField
        (
            interface.phase2(),
            "F",
           -fvc::flux(F),
            Fs
        );
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    {
        const phaseInterface& interface =
            wallLubricationModelIter()->interface();

        const volVectorField F(wallLubricationModelIter()->F());

        addField
        (
            interface.phase1(),
            "F",
            fvc::flux(F),
            Fs
        );
        addField
        (
            interface.phase2(),
            "F",
           -fvc::flux(F),
            Fs
        );
    }

    // Add the phase pressure
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        const phaseModel& phase = fluid_.movingPhases()[movingPhasei];

        addField
        (
            phase,
            "F",
            phase.pPrimef()*fvc::snGrad(phase)*fluid_.mesh().magSf(),
            Fs
        );
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField Df
        (
            fvc::interpolate(turbulentDispersionModelIter()->D())
        );

        const volScalarField alpha12(interface.phase1() + interface.phase2());
        const surfaceScalarField snGradAlpha1By12
        (
            fvc::snGrad
            (
                interface.phase1()
               /max(alpha12, interface.phase1().residualAlpha())
            )*fluid_.mesh().magSf()
        );
        const surfaceScalarField snGradAlpha2By12
        (
            fvc::snGrad
            (
                interface.phase2()
               /max(alpha12, interface.phase2().residualAlpha())
            )*fluid_.mesh().magSf()
        );

        addField(interface.phase1(), "F", Df*snGradAlpha1By12, Fs);
        addField(interface.phase2(), "F", Df*snGradAlpha2By12, Fs);
    }

    return Fs;
}


Foam::PtrList<Foam::surfaceScalarField>
Foam::momentumTransferSystem::Ffs() const
{
    PtrList<surfaceScalarField> Ffs(fluid_.phases().size());

    // Add the lift force
    forAllConstIter
    (
        liftModelTable,
        liftModels_,
        liftModelIter
    )
    {
        const phaseInterface& interface = liftModelIter()->interface();

        const surfaceScalarField Ff(liftModelIter()->Ff());

        addField
        (
            interface.phase1(),
            "Ff",
            Ff,
            Ffs
        );
        addField
        (
            interface.phase2(),
            "Ff",
           -Ff,
            Ffs
        );
    }

    // Add the wall lubrication force
    forAllConstIter
    (
        wallLubricationModelTable,
        wallLubricationModels_,
        wallLubricationModelIter
    )
    {
        const phaseInterface& interface =
            wallLubricationModelIter()->interface();

        const surfaceScalarField Ff(wallLubricationModelIter()->Ff());

        addField
        (
            interface.phase1(),
            "Ff",
            Ff,
            Ffs
        );
        addField
        (
            interface.phase2(),
            "Ff",
           -Ff,
            Ffs
        );
    }

    // Add the phase pressure
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        const phaseModel& phase = fluid_.movingPhases()[movingPhasei];

        addField
        (
            phase,
            "Ff",
            phase.pPrimef()*fvc::snGrad(phase)*fluid_.mesh().magSf(),
            Ffs
        );
    }

    // Add the turbulent dispersion force
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField Df
        (
            fvc::interpolate(turbulentDispersionModelIter()->D())
        );

        const volScalarField alpha12(interface.phase1() + interface.phase2());
        const surfaceScalarField snGradAlpha1By12
        (
            fvc::snGrad
            (
                interface.phase1()
               /max(alpha12, interface.phase1().residualAlpha())
            )*fluid_.mesh().magSf()
        );
        const surfaceScalarField snGradAlpha2By12
        (
            fvc::snGrad
            (
                interface.phase2()
               /max(alpha12, interface.phase2().residualAlpha())
            )*fluid_.mesh().magSf()
        );

        addField(interface.phase1(), "F", Df*snGradAlpha1By12, Ffs);
        addField(interface.phase2(), "F", Df*snGradAlpha2By12, Ffs);
    }

    return Ffs;
}


void Foam::momentumTransferSystem::invADVs
(
    List<UPtrList<scalarField>>& ADVs
) const
{
    const label n = ADVs.size();

    scalarSquareMatrix AD(n);
    scalarField source(n);
    labelList pivotIndices(n);

    forAll(ADVs[0][0], ci)
    {
        for (label i=0; i<n; i++)
        {
            for (label j=0; j<n; j++)
            {
                AD(i, j) = ADVs[i][j][ci];
            }
        }

        // Calculate the inverse of AD using LD decomposition
        // and back-substitution
        LUDecompose(AD, pivotIndices);

        for (label j=0; j<n; j++)
        {
            source = Zero;
            source[j] = 1;

            LUBacksubstitute(AD, pivotIndices, source);

            for (label i=0; i<n; i++)
            {
                ADVs[i][j][ci] = source[i];
            }
        }
    }
}


template<class GeoMesh>
void Foam::momentumTransferSystem::invADVs
(
    PtrList<PtrList<GeometricField<scalar, GeoMesh>>>& ADVs
) const
{
    const label n = ADVs.size();

    List<UPtrList<scalarField>> ADps(n);

    for (label i=0; i<n; i++)
    {
        ADps[i].setSize(n);

        for (label j=0; j<n; j++)
        {
            ADps[i].set(j, &ADVs[i][j]);

            ADVs[i][j].dimensions().reset(dimless/ADVs[i][j].dimensions());
        }
    }

    invADVs(ADps);

    forAll(ADVs[0][0].boundaryField(), patchi)
    {
        for (label i=0; i<n; i++)
        {
            for (label j=0; j<n; j++)
            {
                ADps[i].set(j, &ADVs[i][j].boundaryFieldRef()[patchi]);
            }
        }

        invADVs(ADps);
    }
}


void Foam::momentumTransferSystem::invADVs
(
    const PtrList<volScalarField>& As,
    PtrList<volVectorField>& HVms,
    PtrList<PtrList<volScalarField>>& invADVs,
    PtrList<PtrList<surfaceScalarField>>& invADVfs
) const
{
    const label n = As.size();

    invADVs.setSize(n);
    invADVfs.setSize(n);

    forAll(invADVs, i)
    {
        invADVs.set(i, new PtrList<volScalarField>(n));
        invADVs[i].set(i, As[i].clone());

        invADVfs.set(i, new PtrList<surfaceScalarField>(n));
        invADVfs[i].set(i, fvc::interpolate(As[i]));
    }

    labelList movingPhases(fluid_.phases().size(), -1);
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        movingPhases[fluid_.movingPhases()[movingPhasei].index()] =
            movingPhasei;
    }

    Kds_.clear();

    // Update the drag coefficients
    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        const phaseInterface& interface = dragModelIter()->interface();

        tmp<volScalarField> tKd(dragModelIter()->K());
        volScalarField& Kd = tKd.ref();
        Kds_.insert(dragModelIter.key(), tKd.ptr());

        // Zero-gradient the drag coefficient to boundaries with fixed velocity
        forAll(Kd.boundaryField(), patchi)
        {
            if
            (
                (
                    !interface.phase1().stationary()
                 && interface.phase1().U()()
                   .boundaryField()[patchi].fixesValue()
                )
             && (
                    !interface.phase2().stationary()
                 && interface.phase2().U()()
                   .boundaryField()[patchi].fixesValue()
                )
            )
            {
                Kd.boundaryFieldRef()[patchi] =
                    Kd.boundaryField()[patchi].patchInternalField();
            }
        }

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const volScalarField Kdij
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))*Kd
                );

                const surfaceScalarField Kdijf(fvc::interpolate(Kdij));

                invADVs[i][i] += Kdij;
                invADVfs[i][i] += Kdijf;

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    invADVs[i].set(j, -Kdij);
                    invADVfs[i].set(j, -Kdijf);
                }
            }
        }
    }

    // Clear the Kds_ if they are not needed for the optional dragCorrection
    const pimpleNoLoopControl& pimple = fluid_.pimple();
    if (!pimple.dict().lookupOrDefault<Switch>("dragCorrection", false))
    {
        Kds_.clear();
    }

    for (label i=0; i<n; i++)
    {
        for (label j=0; j<n; j++)
        {
            if (!invADVs[i].set(j))
            {
                invADVs[i].set
                (
                    j,
                    volScalarField::New
                    (
                        "0",
                        fluid_.mesh(),
                        dimensionedScalar(As[0].dimensions(), 0)
                    )
                );

                invADVfs[i].set
                (
                    j,
                    surfaceScalarField::New
                    (
                        "0",
                        fluid_.mesh(),
                        dimensionedScalar(As[0].dimensions(), 0)
                    )
                );
            }
        }
    }

    // Cache the phase acceleration As and Hs
    PtrList<volScalarField> ADUDts(movingPhases.size());
    PtrList<volVectorField> HDUDts(movingPhases.size());

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        VmIter
    )
    {
        const phaseInterface& interface = VmIter()->interface();

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const label i = movingPhases[phase.index()];

            if (i != -1 && !ADUDts.set(i))
            {
                const fvVectorMatrix DUDt(phase.DUDt());
                ADUDts.set(i, DUDt.A());
                HDUDts.set(i, DUDt.H());
            }
        }
    }

    // Add the virtual mass contributions
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        VmIter
    )
    {
        const phaseInterface& interface = VmIter()->interface();
        const volScalarField Vm(VmIter()->K());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))*Vm
                );

                {
                    const volScalarField AVm(VmPhase*ADUDts[i]);

                    invADVs[i][i] += AVm;
                    invADVfs[i][i] += fvc::interpolate(AVm);

                    addField
                    (
                        i,
                        IOobject::groupName("HVm", phase.name()),
                        VmPhase*HDUDts[i],
                        HVms
                    );
                }

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    const volScalarField AVm(VmPhase*ADUDts[j]);

                    invADVs[i][j] -= AVm;
                    invADVfs[i][j] -= fvc::interpolate(AVm);

                    addField
                    (
                        i,
                        IOobject::groupName("HVm", phase.name()),
                        -VmPhase*HDUDts[j],
                        HVms
                    );
                }
            }
        }
    }

    momentumTransferSystem::invADVs(invADVs);
    momentumTransferSystem::invADVs(invADVfs);
}


Foam::PtrList<Foam::PtrList<Foam::surfaceScalarField>>
Foam::momentumTransferSystem::invADVfs
(
    const PtrList<surfaceScalarField>& Afs,
    PtrList<surfaceScalarField>& HVmfs
) const
{
    const label n = Afs.size();

    PtrList<PtrList<surfaceScalarField>> invADVfs(n);

    forAll(invADVfs, i)
    {
        invADVfs.set(i, new PtrList<surfaceScalarField>(n));
        invADVfs[i].set(i, Afs[i].clone());
    }

    labelList movingPhases(fluid_.phases().size(), -1);
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        movingPhases[fluid_.movingPhases()[movingPhasei].index()] =
            movingPhasei;
    }

    forAllConstIter
    (
        dragModelTable,
        dragModels_,
        dragModelIter
    )
    {
        const phaseInterface& interface = dragModelIter()->interface();
        const surfaceScalarField Kdf(dragModelIter()->Kf());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const surfaceScalarField alphaf
                (
                    fvc::interpolate(otherPhase)
                );

                const surfaceScalarField Kdfij
                (
                    (alphaf/max(alphaf, otherPhase.residualAlpha()))*Kdf
                );

                invADVfs[i][i] += Kdfij;

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    invADVfs[i].set(j, -Kdfij);
                }
            }
        }
    }

    for (label i=0; i<n; i++)
    {
        for (label j=0; j<n; j++)
        {
            if (!invADVfs[i].set(j))
            {
                invADVfs[i].set
                (
                    j,
                    surfaceScalarField::New
                    (
                        "0",
                        fluid_.mesh(),
                        dimensionedScalar(Afs[0].dimensions(), 0)
                    )
                );
            }
        }
    }

    // Cache the phase acceleration Afs and Hs
    PtrList<volScalarField> AUgradUs(movingPhases.size());
    PtrList<volVectorField> HUgradUs(movingPhases.size());

    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        VmIter
    )
    {
        const phaseInterface& interface = VmIter()->interface();

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const label i = movingPhases[phase.index()];

            if (i != -1 && !AUgradUs.set(i))
            {
                const fvVectorMatrix UgradU(phase.UgradU());
                AUgradUs.set(i, UgradU.A());
                HUgradUs.set(i, UgradU.H());
            }
        }
    }

    // Add the virtual mass contributions
    forAllConstIter
    (
        virtualMassModelTable,
        virtualMassModels_,
        VmIter
    )
    {
        const phaseInterface& interface = VmIter()->interface();
        const volScalarField Vm(VmIter()->K());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];

            if (i != -1)
            {
                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                   *Vm
                );

                const surfaceScalarField VmPhasef(fvc::interpolate(VmPhase));

                invADVfs[i][i] +=
                    byDt(VmPhasef) + fvc::interpolate(VmPhase*AUgradUs[i]);

                addField
                (
                    i,
                    IOobject::groupName("HVmf", phase.name()),
                    VmPhasef
                   *byDt
                    (
                        fvc::absolute
                        (
                            fluid_.MRF().absolute
                            (
                                phase.phi()().oldTime()
                            ),
                            phase.U()
                        )
                    )
                  + fvc::flux(VmPhase*HUgradUs[i]),
                    HVmfs
                );

                const label j = movingPhases[otherPhase.index()];

                if (j != -1)
                {
                    invADVfs[i][j] -=
                        byDt(VmPhasef) + fvc::interpolate(VmPhase*AUgradUs[j]);

                    addField
                    (
                        i,
                        IOobject::groupName("HVmf", phase.name()),
                       -VmPhasef
                       *byDt
                        (
                            fvc::absolute
                            (
                                fluid_.MRF().absolute
                                (
                                    otherPhase.phi()().oldTime()
                                ),
                                otherPhase.U()
                            )
                        )
                      - fvc::flux(VmPhase*HUgradUs[j]),
                        HVmfs
                    );
                }
            }
        }
    }

    invADVs(invADVfs);

    return invADVfs;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::momentumTransferSystem::alphaDByAf
(
    const PtrList<volScalarField>& rAs
) const
{
    tmp<surfaceScalarField> alphaDByAf;

    // Add the phase pressure
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        const phaseModel& phase = fluid_.movingPhases()[movingPhasei];

        addTmpField
        (
            alphaDByAf,
            fvc::interpolate(max(phase, scalar(0)))
           *fvc::interpolate(rAs[phase.index()])
           *phase.pPrimef()
        );
    }

    // Add the turbulent dispersion
    forAllConstIter
    (
        turbulentDispersionModelTable,
        turbulentDispersionModels_,
        turbulentDispersionModelIter
    )
    {
        const phaseInterface& interface =
            turbulentDispersionModelIter()->interface();

        const surfaceScalarField alpha1f
        (
            fvc::interpolate(max(interface.phase1(), scalar(0)))
        );

        const surfaceScalarField alpha2f
        (
            fvc::interpolate(max(interface.phase2(), scalar(0)))
        );

        addTmpField
        (
            alphaDByAf,
            alpha1f*alpha2f
           /max(alpha1f + alpha2f, interface.phase1().residualAlpha())
           *fvc::interpolate
            (
                max
                (
                    rAs[interface.phase1().index()],
                    rAs[interface.phase2().index()]
                )
               *turbulentDispersionModelIter()->D()
            )
        );
    }

    return alphaDByAf;
}


Foam::PtrList<Foam::surfaceScalarField>
Foam::momentumTransferSystem::ddtCorrs() const
{
    PtrList<surfaceScalarField> ddtCorrs(fluid_.phases().size());

    // Add correction
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        const phaseModel& phase = fluid_.movingPhases()[movingPhasei];

        addField
        (
            phase,
            "ddtCorr",
            fvc::ddtCorr
            (
                phase,
                phase.rho(),
                phase.U()(),
                phase.phi()(),
                phase.Uf()
            ),
            ddtCorrs
        );
    }


    const pimpleNoLoopControl& pimple = fluid_.pimple();
    const Switch VmDdtCorr
    (
        pimple.dict().lookupOrDefault<Switch>("VmDdtCorrection", false)
    );

    // Optionally ddd virtual mass correction
    if (VmDdtCorr)
    {
        PtrList<volScalarField> VmDdtCoeffs(fluid_.phases().size());
        PtrList<surfaceScalarField> VmDdtCorrs(fluid_.phases().size());

        forAll(fluid_.movingPhases(), movingPhasei)
        {
            const phaseModel& phase = fluid_.movingPhases()[movingPhasei];
            const label phasei = phase.index();

            VmDdtCorrs.set
            (
                phasei,
                fvc::ddtCorr
                (
                    phase.U()(),
                    phase.phi()(),
                    phase.Uf()
                )
            );
        }

        forAllConstIter
        (
            virtualMassModelTable,
            virtualMassModels_,
            VmIter
        )
        {
            const phaseInterface& interface = VmIter()->interface();
            const volScalarField Vm(VmIter()->K());

            forAllConstIter(phaseInterface, interface, iter)
            {
                const phaseModel& phase = iter();
                const label phasei = phase.index();
                const phaseModel& otherPhase = iter.otherPhase();
                const label otherPhasei = otherPhase.index();

                const volScalarField VmPhase
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))
                   *Vm
                );

                addField
                (
                    phase,
                    "ddtCorr",
                    fvc::interpolate(VmPhase)
                   *(VmDdtCorrs[phasei] - VmDdtCorrs[otherPhasei]),
                    ddtCorrs
                );
            }
        }
    }

    return ddtCorrs;
}


void Foam::momentumTransferSystem::dragCorrs
(
    PtrList<volVectorField>& dragCorrs,
    PtrList<surfaceScalarField>& dragCorrfs
) const
{
    labelList movingPhases(fluid_.phases().size(), -1);
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        movingPhases[fluid_.movingPhases()[movingPhasei].index()] =
            movingPhasei;
    }

    PtrList<volVectorField> Uphis(fluid_.movingPhases().size());
    forAll(fluid_.movingPhases(), movingPhasei)
    {
        Uphis.set
        (
            movingPhasei,
            fvc::reconstruct(fluid_.movingPhases()[movingPhasei].phi())
        );
    }

    forAllConstIter(KdTable, Kds_, KdIter)
    {
        const volScalarField& K(*KdIter());
        const phaseInterface interface(fluid_, KdIter.key());

        forAllConstIter(phaseInterface, interface, iter)
        {
            const phaseModel& phase = iter();
            const phaseModel& otherPhase = iter.otherPhase();

            const label i = movingPhases[phase.index()];
            const label j = movingPhases[otherPhase.index()];

            if (i != -1)
            {
                const volScalarField K1
                (
                    (otherPhase/max(otherPhase, otherPhase.residualAlpha()))*K
                );

                addField
                (
                    i,
                    IOobject::groupName("dragCorr", phase.name()),
                    K1*(j == -1 ? -Uphis[i] : (Uphis[j] - Uphis[i])),
                    dragCorrs
                );

                addField
                (
                    i,
                    IOobject::groupName("dragCorrf", phase.name()),
                    fvc::interpolate(K1)
                   *(j == -1 ? -phase.phi() : (otherPhase.phi() - phase.phi())),
                    dragCorrfs
                );
            }
        }
    }
}


bool Foam::momentumTransferSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
