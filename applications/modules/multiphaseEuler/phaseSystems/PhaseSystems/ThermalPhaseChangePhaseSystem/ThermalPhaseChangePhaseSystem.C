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

#include "ThermalPhaseChangePhaseSystem.H"
#include "heatTransferModel.H"
#include "alphatPhaseChangeWallFunctionBase.H"
#include "fvcVolumeIntegrate.H"
#include "fvmSup.H"
#include "rhoMulticomponentThermo.H"
#include "wallBoilingHeatTransfer.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::addDmdts
(
    PtrList<volScalarField>& dmdts
) const
{
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        addField(interface.phase1(), "dmdt", *dmdtfIter(), dmdts);
        addField(interface.phase2(), "dmdt", - *dmdtfIter(), dmdts);
    }

    forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
    {
        const phaseInterface interface(*this, nDmdtfIter.key());

        addField(interface.phase1(), "dmdt", *nDmdtfIter(), dmdts);
        addField(interface.phase2(), "dmdt", - *nDmdtfIter(), dmdts);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
ThermalPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    volatile_(this->template lookupOrDefault<word>("volatile", "none")),
    dmdt0s_(this->phases().size()),
    pressureImplicit_
    (
        this->template lookupOrDefault<Switch>("pressureImplicit", true)
    )
{
    this->generateInterfacialModels(saturationModels_);

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phaseInterface& interface = saturationModelIter()->interface();
        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        this->template validateMassTransfer
        <
            interfaceSaturationTemperatureModel
        >(interface);

        if
        (
            !this->heatTransferModels_.found(interface)
         || !this->heatTransferModels_[interface]->haveModelInThe(phase1)
         || !this->heatTransferModels_[interface]->haveModelInThe(phase2)
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for both sides of the "
                 << interface.name() << " interface is not specified. This is "
                 << "required by the corresponding saturation model"
                 << exit(FatalError);
        }
    }

    // Generate interfacial mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phaseInterface& interface = saturationModelIter()->interface();

        dmdtfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:dmdtf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        d2mdtdpfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:d2mdtdpf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar((dimDensity/dimTime)/dimPressure, 0)
            )
        );

        Tfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:Tf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                (
                    interface.phase1().thermo().T()
                  + interface.phase2().thermo().T()
                )/2
            )
        );

        Tsats_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:Tsat",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                saturationModels_[interface]->Tsat
                (
                    interface.phase1().thermo().p()
                )
            )
        );

        nDmdtfs_.insert
        (
            interface,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "thermalPhaseChange:nucleation:dmdtf",
                        interface.name()
                    ),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
~ThermalPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
const Foam::interfaceSaturationTemperatureModel&
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::saturation
(
    const phaseInterfaceKey& key
) const
{
    return saturationModels_[key];
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phaseInterfaceKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    if (dmdtfs_.found(key))
    {
        tDmdtf.ref() += *dmdtfs_[key];
    }

    if (nDmdtfs_.found(key))
    {
        tDmdtf.ref() += *nDmdtfs_[key];
    }

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    addDmdts(dmdts);

    return dmdts;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::d2mdtdps() const
{
    PtrList<volScalarField> d2mdtdps(BasePhaseSystem::d2mdtdps());

    forAllConstIter(phaseSystem::dmdtfTable, d2mdtdpfs_, d2mdtdpfIter)
    {
        const phaseInterface interface(*this, d2mdtdpfIter.key());

        addField(interface.phase1(), "d2mdtdp", *d2mdtdpfIter(), d2mdtdps);
        addField(interface.phase2(), "d2mdtdp", - *d2mdtdpfIter(), d2mdtdps);
    }

    return d2mdtdps;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);
    this->addDmdtUfs(nDmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    // Create temperatures at which to evaluate nucleation mass transfers
    phaseSystem::dmdtfTable Tns;
    forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
    {
        const phaseInterface interface(*this, nDmdtfIter.key());
        const interfaceSaturationTemperatureModel& satModel =
            this->saturation(nDmdtfIter.key());

        Tns.insert
        (
            interface,
            satModel.Tsat(interface.phase1().thermo().p()).ptr()
        );
    }

    // Mass transfer terms
    if (volatile_ != "none")
    {
        {
            phaseSystem::dmidtfTable dmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
            {
                const phaseInterface interface(*this, dmdtfIter.key());

                dmidtfs.insert(interface, new HashPtrTable<volScalarField>());
                dmidtfs[interface]->insert
                (
                    volatile_,
                    new volScalarField(*dmdtfIter())
                );
            }

            this->addDmidtHefs
            (
                dmidtfs,
                //Tfs_,
                Tsats_,
                latentHeatScheme::upwind,
                latentHeatTransfer::mass,
                eqns
            );
        }
        {
            phaseSystem::dmidtfTable nDmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
            {
                const phaseInterface interface(*this, nDmdtfIter.key());

                nDmidtfs.insert(interface, new HashPtrTable<volScalarField>());
                nDmidtfs[interface]->insert
                (
                    volatile_,
                    new volScalarField(*nDmdtfIter())
                );
            }

            this->addDmidtHefs
            (
                nDmidtfs,
                Tns,
                0,
                latentHeatScheme::upwind,
                eqns
            );
        }
    }
    else
    {
        this->addDmdtHefs
        (
            dmdtfs_,
            //Tfs_,
            Tsats_,
            latentHeatScheme::upwind,
            latentHeatTransfer::mass,
            eqns
        );
        this->addDmdtHefs
        (
            nDmdtfs_,
            Tns,
            0,
            latentHeatScheme::upwind,
            eqns
        );
    }

    // Lagging
    {
        PtrList<volScalarField> dmdts(this->phases().size());

        addDmdts(dmdts);

        forAll(this->phases(), phasei)
        {
            const phaseModel& phase = this->phases()[phasei];

            if (dmdt0s_.set(phase.index()))
            {
                *eqns[phase.name()] +=
                    fvm::Sp
                    (
                        dmdt0s_[phase.index()] - dmdts[phase.index()],
                        phase.thermo().he()
                    );
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    if (volatile_ != "none")
    {
        {
            phaseSystem::dmidtfTable dmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
            {
                const phaseInterface interface(*this, dmdtfIter.key());

                dmidtfs.insert(interface, new HashPtrTable<volScalarField>());
                dmidtfs[interface]->insert
                (
                    volatile_,
                    new volScalarField(*dmdtfIter())
                );
            }

            this->addDmidtYf(dmidtfs, eqns);
        }

        {
            phaseSystem::dmidtfTable nDmidtfs;

            forAllConstIter(phaseSystem::dmdtfTable, nDmdtfs_, nDmdtfIter)
            {
                const phaseInterface interface(*this, nDmdtfIter.key());

                nDmidtfs.insert(interface, new HashPtrTable<volScalarField>());
                nDmidtfs[interface]->insert
                (
                    volatile_,
                    new volScalarField(*nDmdtfIter())
                );
            }

            this->addDmidtYf(nDmidtfs, eqns);
        }
    }
    else
    {
        this->addDmdtYfs(dmdtfs_, eqns);
        this->addDmdtYfs(nDmdtfs_, eqns);
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::
correctContinuityError()
{
    dmdt0s_ = PtrList<volScalarField>(this->phases().size());

    addDmdts(dmdt0s_);

    BasePhaseSystem::correctContinuityError();
}


template<class BasePhaseSystem>
void
Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::correctInterfaceThermo()
{
    typedef
        Foam::heatTransferModels::wallBoilingHeatTransfer
        wallBoilingHeatTransferModel;

    HashTable<const wallBoilingHeatTransferModel*>
        wallBoilingHeatTransferModels =
        this->mesh().template lookupClass<wallBoilingHeatTransferModel>();

    typedef
        compressible::alphatPhaseChangeWallFunctionBase
        alphatPhaseChangeWallFunction;

    forAllConstIter
    (
        saturationModelTable,
        saturationModels_,
        saturationModelIter
    )
    {
        const phaseInterface& interface = saturationModelIter()->interface();
        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();
        const volScalarField& T1(thermo1.T());
        const volScalarField& T2(thermo2.T());

        const sidedBlendedHeatTransferModel& heatTransferModel =
            this->heatTransferModels_[interface];

        // Interfacial mass transfer update
        {
            volScalarField& dmdtf(*this->dmdtfs_[interface]);
            volScalarField& Tf(*this->Tfs_[interface]);

            volScalarField& Tsat(*this->Tsats_[interface]);
            Tsat = saturationModelIter()->Tsat(thermo1.p());

            const volScalarField L
            (
                volatile_ != "none"
              ? this->Li
                (
                    interface,
                    volatile_,
                    dmdtf,
                    Tsat,
                    latentHeatScheme::symmetric
                )
              : this->L
                (
                    interface,
                    dmdtf,
                    Tsat,
                    latentHeatScheme::symmetric
                )
            );

            volScalarField H1(heatTransferModel.modelInThe(phase1).K(0));
            volScalarField H2(heatTransferModel.modelInThe(phase2).K(0));

            volScalarField dmdtfNew((H1*(Tsat - T1) + H2*(Tsat - T2))/L);

            if (volatile_ != "none")
            {
                dmdtfNew *=
                    neg0(dmdtfNew)*phase1.Y(volatile_)
                  + pos(dmdtfNew)*phase2.Y(volatile_);
            }

            if (pressureImplicit_)
            {
                volScalarField& d2mdtdpf(*this->d2mdtdpfs_[interface]);

                const dimensionedScalar dp(rootSmall*thermo1.p().average());

                const volScalarField dTsatdp
                (
                    (
                        saturationModelIter()->Tsat(thermo1.p() + dp/2)
                      - saturationModelIter()->Tsat(thermo1.p() - dp/2)
                    )/dp
                );

                d2mdtdpf = (H1 + H2)*dTsatdp/L;

                if (volatile_ != "none")
                {
                    d2mdtdpf *=
                        neg0(dmdtfNew)*phase1.Y(volatile_)
                      + pos(dmdtfNew)*phase2.Y(volatile_);
                }
            }

            H1 = heatTransferModel.modelInThe(phase1).K();
            H2 = heatTransferModel.modelInThe(phase2).K();

            // Limit the H[12] to avoid /0
            H1.max(small);
            H2.max(small);

            Tf = (H1*T1 + H2*T2 + dmdtfNew*L)/(H1 + H2);

            Info<< Tsat.name()
                << ": min = " << gMin(Tsat.primitiveField())
                << ", mean = " << gAverage(Tsat.primitiveField())
                << ", max = " << gMax(Tsat.primitiveField())
                << endl;

            Info<< Tf.name()
                << ": min = " << gMin(Tf.primitiveField())
                << ", mean = " << gAverage(Tf.primitiveField())
                << ", max = " << gMax(Tf.primitiveField())
                << endl;

            const scalar dmdtfRelax =
                this->mesh().solution().fieldRelaxationFactor(dmdtf.member());

            dmdtf = (1 - dmdtfRelax)*dmdtf + dmdtfRelax*dmdtfNew;

            Info<< dmdtf.name()
                << ": min = " << gMin(dmdtf.primitiveField())
                << ", mean = " << gAverage(dmdtf.primitiveField())
                << ", max = " << gMax(dmdtf.primitiveField())
                << ", integral = " << fvc::domainIntegrate(dmdtf).value()
                << endl;
        }

        // Nucleation mass transfer update
        {
            volScalarField& nDmdtf(*this->nDmdtfs_[interface]);
            nDmdtf = Zero;

            bool wallBoilingActive = false;

            forAllConstIter
            (
                HashTable<const wallBoilingHeatTransferModel*>,
                wallBoilingHeatTransferModels,
                wallBoilingHeatTransferModelIter
            )
            {
                const wallBoilingHeatTransferModel& wbht =
                    *wallBoilingHeatTransferModelIter();

                if (!wbht.activePhaseInterface(interface)) continue;

                wallBoilingActive = true;

                nDmdtf +=
                    (interface == wbht.activePhaseInterface() ? +1 : -1)
                   *wbht.dmdtf();
            }

            forAllConstIter(phaseInterface, interface, interfaceIter)
            {
                const phaseModel& phase = interfaceIter();

                const word alphatName =
                    IOobject::groupName("alphat", phase.name());

                if (!phase.mesh().foundObject<volScalarField>(alphatName))
                    continue;

                const volScalarField& alphat =
                    phase.mesh().lookupObject<volScalarField>(alphatName);

                forAll(alphat.boundaryField(), patchi)
                {
                    const fvPatchScalarField& alphatp =
                        alphat.boundaryField()[patchi];

                    if (!isA<alphatPhaseChangeWallFunction>(alphatp)) continue;

                    const alphatPhaseChangeWallFunction& alphatw =
                        refCast<const alphatPhaseChangeWallFunction>(alphatp);

                    if (!alphatw.activeInterface(interface)) continue;

                    wallBoilingActive = true;

                    UIndirectList<scalar> nDmdtfp
                    (
                        nDmdtf.primitiveFieldRef(),
                        alphatp.patch().faceCells()
                    );

                    nDmdtfp =
                        scalarField(nDmdtfp)
                      - (interfaceIter.index() == 0 ? +1 : -1)
                       *alphatw.dmdtf();
                }
            }

            if (wallBoilingActive)
            {
                Info<< nDmdtf.name()
                    << ": min = " << gMin(nDmdtf.primitiveField())
                    << ", mean = " << gAverage(nDmdtf.primitiveField())
                    << ", max = " << gMax(nDmdtf.primitiveField())
                    << ", integral = " << fvc::domainIntegrate(nDmdtf).value()
                    << endl;
            }
        }
    }
}


template<class BasePhaseSystem>
bool Foam::ThermalPhaseChangePhaseSystem<BasePhaseSystem>::read()
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
