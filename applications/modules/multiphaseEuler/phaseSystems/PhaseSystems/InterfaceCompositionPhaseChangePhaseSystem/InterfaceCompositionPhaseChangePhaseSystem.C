/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "interfaceCompositionModel.H"
#include "heatTransferModel.H"
#include "diffusiveMassTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctDmdtfs()
{
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = compositionModel.interface();

        *dmdtfs_[interface] = Zero;

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            const phaseModel& phase = interfaceIter();

            if (!compositionModel.haveModelInThe(phase)) continue;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.modelInThe(phase).species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                *dmdtfs_[interface] +=
                    (interfaceIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[interface])[specie]
                      + *(*dmidtfSps_[interface])[specie]*phase.Y(specie)
                    );
            }
        }
    }
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmidtfTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
totalDmidtfs() const
{
    autoPtr<phaseSystem::dmidtfTable> totalDmidtfsPtr
    (
        new phaseSystem::dmidtfTable
    );
    phaseSystem::dmidtfTable& totalDmidtfs = totalDmidtfsPtr();

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = compositionModel.interface();

        if (!totalDmidtfs.found(interface))
        {
            totalDmidtfs.insert(interface, new HashPtrTable<volScalarField>());
        }

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            const phaseModel& phase = interfaceIter();

            if (!compositionModel.haveModelInThe(phase)) continue;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.modelInThe(phase).species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                tmp<volScalarField> dmidtf
                (
                    (interfaceIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[interface])[specie]
                      + *(*dmidtfSps_[interface])[specie]*phase.Y(specie)
                    )
                );

                if (totalDmidtfs[interface]->found(specie))
                {
                    *(*totalDmidtfs[interface])[specie] += dmidtf;
                }
                else
                {
                    totalDmidtfs[interface]->insert(specie, dmidtf.ptr());
                }
            }
        }
    }

    return totalDmidtfsPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
InterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    nInterfaceCorrectors_
    (
        this->template lookupOrDefault<label>("nInterfaceCorrectors", 1)
    )
{
    this->generateInterfacialModels(interfaceCompositionModels_);
    this->generateInterfacialModels(diffusiveMassTransferModels_);

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& sidedCompositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = sidedCompositionModel.interface();
        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        this->template
            validateMassTransfer<interfaceCompositionModel>(interface);

        if (!this->diffusiveMassTransferModels_.found(interface))
        {
            FatalErrorInFunction
                << "A diffusive mass transfer model for the " << interface
                << " interface is not specified. This is required by the "
                << "corresponding interface composition model."
                << exit(FatalError);
        }

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            if
            (
                sidedCompositionModel.haveModelInThe(interfaceIter())
             && !diffusiveMassTransferModels_[interface]
                ->haveModelInThe(interfaceIter())
            )
            {
                FatalErrorInFunction
                    << "A diffusive mass transfer model for the "
                    << interfaceIter().name() << " side of the "
                    << interface.name() << " interface is not "
                    << "specified. This is required by the corresponding "
                    << "interface composition model."
                    << exit(FatalError);
            }
        }

        if
        (
            !this->heatTransferModels_.found(interface)
         || !this->heatTransferModels_[interface]->haveModelInThe(phase1)
         || !this->heatTransferModels_[interface]->haveModelInThe(phase2)
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for both sides of the "
                 << interface.name()
                 << " interface is not specified. This is required by the "
                 << "corresponding interface composition model"
                 << exit(FatalError);
        }
    }

    // Generate mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& sidedCompositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = sidedCompositionModel.interface();

        dmdtfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "interfaceCompositionPhaseChange:dmdtf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, Zero)
            )
        );

        dmidtfSus_.insert(interface, new HashPtrTable<volScalarField>());

        dmidtfSps_.insert(interface, new HashPtrTable<volScalarField>());

        Tfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "interfaceCompositionPhaseChange:Tf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (
                    interface.phase1().thermo().T()
                  + interface.phase2().thermo().T()
                )/2
            )
        );

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            const phaseModel& phase = interfaceIter();

            if (!sidedCompositionModel.haveModelInThe(phase)) continue;

            const interfaceCompositionModel& compositionModel =
                sidedCompositionModel.modelInThe(phase);

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                dmidtfSus_[interface]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSu",
                                    specie
                                ),
                                interface.name()
                            ),
                            this->mesh().time().name(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );

                dmidtfSps_[interface]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSp",
                                    specie
                                ),
                                interface.name()
                            ),
                            this->mesh().time().name(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phaseInterfaceKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    if (dmdtfs_.found(key))
    {
        tDmdtf.ref() += *dmdtfs_[key];
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        addField(interface.phase1(), "dmdt", *dmdtfIter(), dmdts);
        addField(interface.phase2(), "dmdt", - *dmdtfIter(), dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmidtHefs
    (
        totalDmidtfs(),
        Tfs_,
        latentHeatScheme::symmetric,
        latentHeatTransfer::mass,
        eqns
    );

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    // Explicit
    /*
    this->addDmidtYf(totalDmidtfs(), eqns);
    */

    // Semi-implicit
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = compositionModel.interface();

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            const phaseModel& phase = interfaceIter();
            const phaseModel& otherPhase = interfaceIter.otherPhase();

            if (!compositionModel.haveModelInThe(phase)) continue;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.modelInThe(phase).species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                // Implicit transport through this phase
                *eqns[phase.Y(specie).name()] +=
                    *(*dmidtfSus_[interface])[specie]
                  + fvm::Sp(*(*dmidtfSps_[interface])[specie], phase.Y(specie));

                // Explicit transport out of the other phase
                if (eqns.found(IOobject::groupName(specie, otherPhase.name())))
                {
                    *eqns[otherPhase.Y(specie).name()] -=
                        *(*dmidtfSus_[interface])[specie]
                      + *(*dmidtfSps_[interface])[specie]*phase.Y(specie);
                }
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();

    // Sum up the contribution from each interface composition model
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const sidedInterfaceCompositionModel& compositionModel =
            interfaceCompositionModelIter();

        const phaseInterface& interface = compositionModel.interface();

        const volScalarField& Tf(*this->Tfs_[interface]);

        forAllConstIter(phaseInterface, interface, interfaceIter)
        {
            const phaseModel& phase = interfaceIter();

            if (!compositionModel.haveModelInThe(phase)) continue;

            const volScalarField K
            (
                diffusiveMassTransferModels_[interface]->modelInThe(phase).K()
            );

            forAllConstIter
            (
                hashedWordList,
                compositionModel.modelInThe(phase).species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                const volScalarField KD
                (
                    K*compositionModel.modelInThe(phase).D(specie)
                );
                const volScalarField Yf
                (
                    compositionModel.modelInThe(phase).Yf(specie, Tf)
                );

                *(*dmidtfSus_[interface])[specie] = phase.rho()*KD*Yf;
                *(*dmidtfSps_[interface])[specie] = - phase.rho()*KD;
            }
        }
    }

    correctDmdtfs();
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctSpecies()
{
    BasePhaseSystem::correctSpecies();

    correctDmdtfs();
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    // This loop solves for the interface temperatures, Tf, and updates the
    // interface composition models.
    //
    // In the presence of thermally-coupled mass transfer, the relation between
    // heat transfers across the interface between phases 1 and 2 is:
    //
    //                         Q1 + Q2 = mDot*L
    //     H1*(Tf - T1) + H2*(Tf - T1) = K*rho*(Yfi - Yi)*Li
    //
    // Where Q1 and Q2 are the net transfer into phases 1 and 2 respectively,
    // H1 and H2 are the heat transfer coefficients on either side, Tf is the
    // temperature at the interface, mDot is the mass transfer rate from phase
    // 2 to phase 1, and L is the latent heat of phase 2 minus phase 1, K is
    // the diffusive mass transfer coefficient, Yfi - Yi is the concentration
    // difference of a transferring specie between the interface and the bulk
    // driving the transfer, Li is the latent heat change of the specie, and
    // rho is the density in the phase in which the diffusive mass transfer is
    // being represented.
    //
    // Yfi is likely to be a strong non-linear (typically exponential) function
    // of Tf, so the solution for the temperature is newton-accelerated.

    forAllIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        sidedInterfaceCompositionModel& compositionModel =
            *interfaceCompositionModelIter();

        const phaseInterface& interface = compositionModel.interface();
        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        const sidedBlendedHeatTransferModel& heatTransferModel =
            this->heatTransferModels_[interface];

        const sidedBlendedDiffusiveMassTransferModel&
            diffusiveMassTransferModel =
            diffusiveMassTransferModels_[interface];

        const volScalarField H1(heatTransferModel.modelInThe(phase1).K());
        const volScalarField H2(heatTransferModel.modelInThe(phase2).K());
        const dimensionedScalar HSmall("small", heatTransferModel::dimK, small);

        volScalarField& Tf = *this->Tfs_[interface];

        for (label i = 0; i < nInterfaceCorrectors_; ++ i)
        {
            tmp<volScalarField> dmdtLf =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLf", interface.name()),
                    this->mesh(),
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
                );
            tmp<volScalarField> dmdtLfPrime =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLfPrime", interface.name()),
                    this->mesh(),
                    dimensionedScalar(dmdtLf().dimensions()/dimTemperature, 0)
                );

            // Add latent heats from forward and backward models
            forAllConstIter(phaseInterface, interface, interfaceIter)
            {
                const phaseModel& phase = interfaceIter();

                if (!compositionModel.haveModelInThe(phase)) continue;

                const label sign = interfaceIter.index() == 0 ? 1 : -1;

                forAllConstIter
                (
                    hashedWordList,
                    compositionModel.modelInThe(phase).species(),
                    specieIter
                )
                {
                    const word& specie = *specieIter;

                    const volScalarField dY
                    (
                        compositionModel.modelInThe(phase).dY(specie, Tf)
                    );

                    const volScalarField dYfPrime
                    (
                        compositionModel.modelInThe(phase).dYfPrime(specie, Tf)
                    );

                    const volScalarField rhoKDL
                    (
                        phase.rho()
                       *diffusiveMassTransferModel.modelInThe(phase).K()
                       *compositionModel.modelInThe(phase).D(specie)
                       *this->Li
                        (
                            interface,
                            specie,
                            dY,
                            Tf,
                            latentHeatScheme::symmetric
                        )
                    );

                    dmdtLf.ref() += sign*rhoKDL*dY;
                    dmdtLfPrime.ref() += sign*rhoKDL*dYfPrime;
                }
            }

            // Update the interface temperature by applying one step of newton's
            // method to the interface relation
            Tf -=
                (
                    H1*(Tf - interface.phase1().thermo().T())
                  + H2*(Tf - interface.phase2().thermo().T())
                  - dmdtLf
                )
               /(
                    max(H1 + H2 - dmdtLfPrime, HSmall)
                );

            Tf.correctBoundaryConditions();

            Info<< "Tf." << interface.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;

            // Update the interface compositions
            forAllConstIter(phaseInterface, interface, interfaceIter)
            {
                const phaseModel& phase = interfaceIter();

                if (!compositionModel.haveModelInThe(phase)) continue;

                compositionModel.modelInThe(phase).update(Tf);
            }
        }
    }
}


template<class BasePhaseSystem>
bool Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
