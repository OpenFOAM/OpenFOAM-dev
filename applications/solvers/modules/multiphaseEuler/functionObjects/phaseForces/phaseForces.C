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
    phase_
    (
        mesh_.lookupObject<phaseModel>
        (
            IOobject::groupName("alpha", dict.lookup("phase"))
        )
    ),
    fluid_(mesh_.lookupObject<phaseSystem>(phaseSystem::propertiesName))
{
    read(dict);

    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];

        if (&otherPhase == &phase_) continue;

        const phaseInterface interface(phase_, otherPhase);

        if (fluid_.foundInterfacialModel<blendedDragModel>(interface))
        {
            forceFields_.insert
            (
                dragModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName("dragForce", phase_.name()),
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedVector(dimForce/dimVolume, Zero)
                )
            );
        }

        if (fluid_.foundInterfacialModel<blendedVirtualMassModel>(interface))
        {
            forceFields_.insert
            (
                virtualMassModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "virtualMassForce",
                            phase_.name()
                        ),
                        mesh_.time().name(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedVector(dimForce/dimVolume, Zero)
                )
            );
        }

        if (fluid_.foundInterfacialModel<blendedLiftModel>(interface))
        {
            forceFields_.insert
            (
                liftModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName("liftForce", phase_.name()),
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
            fluid_.foundInterfacialModel
            <blendedWallLubricationModel>(interface)
        )
        {
            forceFields_.insert
            (
                wallLubricationModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "wallLubricationForce",
                            phase_.name()
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
            fluid_.foundInterfacialModel
            <blendedTurbulentDispersionModel>(interface)
        )
        {
            forceFields_.insert
            (
                turbulentDispersionModel::typeName,
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "turbulentDispersionForce",
                            phase_.name()
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
    // Zero the force fields
    forAllConstIter
    (
        HashPtrTable<volVectorField>,
        forceFields_,
        forceFieldIter
    )
    {
        *forceFieldIter() = Zero;
    }

    // Add the forces from all the interfaces which contain this phase
    forAll(fluid_.phases(), phasei)
    {
        const phaseModel& otherPhase = fluid_.phases()[phasei];

        if (&otherPhase == &phase_) continue;

        const phaseInterface interface(phase_, otherPhase);

        if (fluid_.foundInterfacialModel<blendedDragModel>(interface))
        {
            *forceFields_[dragModel::typeName] +=
                fluid_.lookupInterfacialModel<blendedDragModel>(interface).K()
               *(otherPhase.U() - phase_.U());
        }

        if (fluid_.foundInterfacialModel<blendedVirtualMassModel>(interface))
        {
            *forceFields_[virtualMassModel::typeName] +=
                fluid_.lookupInterfacialModel
                <blendedVirtualMassModel>(interface).K()
               *(otherPhase.DUDt() - phase_.DUDt());
        }

        if (fluid_.foundInterfacialModel<blendedLiftModel>(interface))
        {
            *forceFields_[liftModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
               *fluid_.lookupInterfacialModel<blendedLiftModel>(interface).F();
        }

        if
        (
            fluid_.foundInterfacialModel
            <blendedWallLubricationModel>(interface)
        )
        {
            *forceFields_[wallLubricationModel::typeName] +=
                (&interface.phase1() == &phase_ ? -1 : +1)
               *fluid_.lookupInterfacialModel
                <blendedWallLubricationModel>(interface).F();
        }

        if
        (
            fluid_.foundInterfacialModel
            <blendedTurbulentDispersionModel>(interface)
        )
        {
            *forceFields_[turbulentDispersionModel::typeName] +=
                fluid_.lookupInterfacialModel
                <blendedTurbulentDispersionModel>(interface).D()
               *fvc::grad
                (
                    otherPhase
                   /max(phase_ + otherPhase, otherPhase.residualAlpha())
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
        forceFields_,
        forceFieldIter
    )
    {
        writeObject(forceFieldIter()->name());
    }

    return true;
}


// ************************************************************************* //
