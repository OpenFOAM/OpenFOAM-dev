/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * * //

template<class modelType>
Foam::tmp<Foam::volVectorField>
Foam::functionObjects::phaseForces::nonDragForce(const phasePair& pair) const
{
    const BlendedInterfacialModel<modelType>& model =
        fluid_.lookupBlendedSubModel<modelType>(pair);

    if (&pair.phase1() == &phase_)
    {
        return -model.template F<vector>();
    }
    else
    {
        return model.template F<vector>();
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
    fluid_(mesh_.lookupObject<phaseSystem>("phaseProperties"))
{
    read(dict);

    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        pairIter
    )
    {
        const phasePair& pair = pairIter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            if (fluid_.foundBlendedSubModel<dragModel>(pair))
            {
                forceFields_.insert
                (
                    dragModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("dragForce", phase_.name()),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid_.foundBlendedSubModel<virtualMassModel>(pair))
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
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid_.foundBlendedSubModel<liftModel>(pair))
            {
                forceFields_.insert
                (
                    liftModel::typeName,
                    new volVectorField
                    (
                        IOobject
                        (
                            IOobject::groupName("liftForce", phase_.name()),
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid_.foundBlendedSubModel<wallLubricationModel>(pair))
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
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }

            if (fluid_.foundBlendedSubModel<turbulentDispersionModel>(pair))
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
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_,
                        dimensionedVector(dimForce/dimVolume, Zero)
                    )
                );
            }
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
    forAllConstIter
    (
        phaseSystem::phasePairTable,
        fluid_.phasePairs(),
        pairIter
    )
    {
        const phasePair& pair = pairIter();

        if (pair.contains(phase_) && !pair.ordered())
        {
            if (fluid_.foundBlendedSubModel<dragModel>(pair))
            {
                *forceFields_[dragModel::typeName] +=
                    fluid_.lookupBlendedSubModel<dragModel>(pair).K()
                   *(pair.otherPhase(phase_).U() - phase_.U());
            }

            if (fluid_.foundBlendedSubModel<virtualMassModel>(pair))
            {
                *forceFields_[virtualMassModel::typeName] +=
                    fluid_.lookupBlendedSubModel<virtualMassModel>(pair).K()
                   *(pair.otherPhase(phase_).DUDt() - phase_.DUDt());
            }

            if (fluid_.foundBlendedSubModel<liftModel>(pair))
            {
                *forceFields_[liftModel::typeName] +=
                    nonDragForce<liftModel>(pair);
            }

            if (fluid_.foundBlendedSubModel<wallLubricationModel>(pair))
            {
                *forceFields_[wallLubricationModel::typeName] +=
                    nonDragForce<wallLubricationModel>(pair);
            }

            if (fluid_.foundBlendedSubModel<turbulentDispersionModel>(pair))
            {
                *forceFields_[turbulentDispersionModel::typeName] +=
                    nonDragForce<turbulentDispersionModel>(pair);
            }
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
