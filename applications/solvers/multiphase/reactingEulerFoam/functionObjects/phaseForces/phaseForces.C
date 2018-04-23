/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
#include "BlendedInterfacialModel.H"
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

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::dragForce() const
{
    tmp<volVectorField> tdragForce
    (
        new volVectorField
        (
            IOobject
            (
                "dragForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& dragForce = tdragForce.ref();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        iter
    )
    {
        const phasePair& pair = iter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            const BlendedInterfacialModel<dragModel>& drag =
                fluid_.lookupBlendedSubModel<dragModel>(pair);

            dragForce += drag.K()*(pair.otherPhase(phase_).U() - phase_.U());
        }
    }

    return tdragForce;
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::virtualMassForce() const
{
    tmp<volVectorField> tvirtualMassForce
    (
        new volVectorField
        (
            IOobject
            (
                "virtualMassForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& virtualMassForce = tvirtualMassForce.ref();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        iter
    )
    {
        const phasePair& pair = iter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            const BlendedInterfacialModel<virtualMassModel>& virtualMass =
                fluid_.lookupBlendedSubModel<virtualMassModel>(pair);

            virtualMassForce +=
                virtualMass.K()
               *(
                    pair.otherPhase(phase_).DUDt()
                  - phase_.DUDt()
                );
        }
    }

    return tvirtualMassForce;
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::liftForce() const
{
    tmp<volVectorField> tLiftForce
    (
        new volVectorField
        (
            IOobject
            (
                "liftForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& liftForce = tLiftForce.ref();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        iter
    )
    {
        const phasePair& pair = iter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            const BlendedInterfacialModel<liftModel>& lift =
                fluid_.lookupBlendedSubModel<liftModel>(pair);

            if (&pair.phase1() == &phase_)
            {
                liftForce += lift.F<vector>();
            }
            else
            {
                liftForce -= lift.F<vector>();
            }
        }
    }

    return tLiftForce;
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::wallLubricationForce() const
{
    tmp<volVectorField> tWallLubricationForce
    (
        new volVectorField
        (
            IOobject
            (
                "wallLubricationForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& wallLubricationForce = tWallLubricationForce.ref();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        iter
    )
    {
        const phasePair& pair = iter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            const BlendedInterfacialModel<wallLubricationModel>&
                wallLubrication =
                    fluid_.lookupBlendedSubModel<wallLubricationModel>(pair);

            if (&pair.phase1() == &phase_)
            {
                wallLubricationForce += wallLubrication.F<vector>();
            }
            else
            {
                wallLubricationForce -= wallLubrication.F<vector>();
            }
        }
    }

    return tWallLubricationForce;
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::turbulentDispersionForce() const
{
    tmp<volVectorField> tturbulentDispersionForce
    (
        new volVectorField
        (
            IOobject
            (
                "turbulentDispersionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimForce/dimVolume, Zero)
        )
    );

    volVectorField& turbulentDispersionForce = tturbulentDispersionForce.ref();

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        iter
    )
    {
        const phasePair& pair = iter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            const BlendedInterfacialModel<turbulentDispersionModel>&
                turbulentDispersion = fluid_.lookupBlendedSubModel
                    <
                        turbulentDispersionModel
                    >(pair);

            if (&pair.phase1() == &phase_)
            {
                turbulentDispersionForce += turbulentDispersion.F<vector>();
            }
            else
            {
                turbulentDispersionForce -= turbulentDispersion.F<vector>();
            }
        }
    }

    return tturbulentDispersionForce;
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
            IOobject::groupName("alpha", dict.lookup("phaseName"))
        )
    ),
    fluid_(mesh_.lookupObject<phaseSystem>("phaseProperties")),
    dragForce_
    (
        IOobject
        (
            IOobject::groupName("dragForce", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimForce/dimVolume, Zero)
    ),
    virtualMassForce_
    (
        IOobject
        (
            IOobject::groupName("virtualMassForce", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimForce/dimVolume, Zero)
    ),
    liftForce_
    (
        IOobject
        (
            IOobject::groupName("liftForce", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimForce/dimVolume, Zero)
    ),
    wallLubricationForce_
    (
        IOobject
        (
            IOobject::groupName("wallLubricationForce", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimForce/dimVolume, Zero)
    ),
    turbulentDispersionForce_
    (
        IOobject
        (
            IOobject::groupName("turbulentDispersionForce", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimForce/dimVolume, Zero)
    )
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
    dragForce_ = dragForce();
    virtualMassForce_ = virtualMassForce();
    liftForce_ = liftForce();
    wallLubricationForce_ = wallLubricationForce();
    turbulentDispersionForce_ = turbulentDispersionForce();

    return true;
}


bool Foam::functionObjects::phaseForces::write()
{
    dragForce_.write();
    virtualMassForce_.write();
    liftForce_.write();
    wallLubricationForce_.write();
    turbulentDispersionForce_.write();

    return true;
}


// ************************************************************************* //
