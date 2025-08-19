/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "phaseForces.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phaseForces, 0);
    addToRunTimeSelectionTable(functionObject, phaseForces, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phaseForces::phaseForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookup<word>("phase")),
    forceFields_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::phaseForces::~phaseForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phaseForces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::phaseForces::execute()
{
    const phaseSystem& fluid =
        mesh_.lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const phaseModel& phase = fluid.phases()[phaseName_];

    // Construct the force fields if needed
    if (!forceFields_.valid())
    {
        forceFields_.set(new HashPtrTable<volVectorField>());

        forAll(fluid.phases(), phasei)
        {
            const phaseModel& otherPhase = fluid.phases()[phasei];

            if (&otherPhase == &phase) continue;

            const phaseInterface interface(phase, otherPhase);

            if (fluid.foundInterfacialModel<blendedDragModel>(interface))
            {
                forceFields_().insert
                (
                    dragModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("dragForce", phase.name()),
                            mesh_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid.foundInterfacialModel<blendedVirtualMassModel>(interface))
            {
                forceFields_().insert
                (
                    virtualMassModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "virtualMassForce",
                                phase.name()
                            ),
                            mesh_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid.foundInterfacialModel<blendedLiftModel>(interface))
            {
                forceFields_().insert
                (
                    liftModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("liftForce", phase.name()),
                            mesh_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                fluid.foundInterfacialModel
                <
                    blendedWallLubricationModel
                >(interface)
            )
            {
                forceFields_().insert
                (
                    wallLubricationModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "wallLubricationForce",
                                phase.name()
                            ),
                            mesh_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if
            (
                fluid.foundInterfacialModel
                <
                    blendedTurbulentDispersionModel
                >(interface)
            )
            {
                forceFields_().insert
                (
                    turbulentDispersionModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "turbulentDispersionForce",
                                phase.name()
                            ),
                            mesh_.time().name(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }
        }
    }

    // Zero the force fields
    forAllConstIter
    (
        HashPtrTable<volVectorField>,
        forceFields_(),
        forceFieldIter
    )
    {
        *forceFieldIter() = Zero;
    }

    // Add the forces from all the interfaces which contain this phase
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid.phases()[phasei];

        if (&otherPhase == &phase) continue;

        const phaseInterface interface(phase, otherPhase);

        if (fluid.foundInterfacialModel<blendedDragModel>(interface))
        {
            *forceFields_()[dragModel::typeName] +=
                fluid.lookupInterfacialModel<blendedDragModel>(interface).K()
               *(otherPhase.U() - phase.U());
        }

        if (fluid.foundInterfacialModel<blendedVirtualMassModel>(interface))
        {
            *forceFields_()[virtualMassModel::typeName] +=
                fluid.lookupInterfacialModel
                <
                    blendedVirtualMassModel
                >(interface).K()
               *(
                    (otherPhase.DUDt() & otherPhase.U())
                  - (phase.DUDt() & phase.U())
                );
        }

        if (fluid.foundInterfacialModel<blendedLiftModel>(interface))
        {
            *forceFields_()[liftModel::typeName] +=
                (&interface.phase1() == &phase ? -1 : +1)
               *fluid.lookupInterfacialModel<blendedLiftModel>(interface).F();
        }

        if
        (
            fluid.foundInterfacialModel
            <
                blendedWallLubricationModel
            >(interface)
        )
        {
            *forceFields_()[wallLubricationModel::typeName] +=
                (&interface.phase1() == &phase ? -1 : +1)
               *fluid.lookupInterfacialModel
                <
                    blendedWallLubricationModel
                >(interface).F();
        }

        if
        (
            fluid.foundInterfacialModel
            <
                blendedTurbulentDispersionModel
            >(interface)
        )
        {
            *forceFields_()[turbulentDispersionModel::typeName] +=
                fluid.lookupInterfacialModel
                <
                    blendedTurbulentDispersionModel
                >(interface).D()
               *fvc::grad
                (
                    otherPhase
                   /max(phase + otherPhase, otherPhase.residualAlpha())
                );
        }
    }

    return true;
}


bool Foam::functionObjects::phaseForces::write()
{
    forAllConstIter
    (
        HashPtrTable<volVectorField>,
        forceFields_(),
        forceFieldIter
    )
    {
        writeObject(forceFieldIter()->name());
    }

    return true;
}


// ************************************************************************* //
