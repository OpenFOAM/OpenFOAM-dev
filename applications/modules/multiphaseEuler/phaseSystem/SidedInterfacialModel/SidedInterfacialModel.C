/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "SidedInterfacialModel.H"
#include "phaseSystem.H"
#include "sidedPhaseInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ModelType>
Foam::SidedInterfacialModel<ModelType>::SidedInterfacialModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.mesh().time().name(),
            interface.mesh()
        )
    ),
    interface_(interface)
{
    // Construct the models
    PtrList<phaseInterface> interfaces;
    PtrList<ModelType> models;
    interface.fluid().generateInterfacialModels<ModelType, sidedPhaseInterface>
    (
        dict,
        interface,
        interfaces,
        models
    );

    // Unpack the interface and model lists to populate the models used for
    // either side of the interface
    forAll(interfaces, i)
    {
        const sidedPhaseInterface& interface =
            refCast<const sidedPhaseInterface>(interfaces[i]);

        if (interface_.index(interface.phase()) == 0)
        {
            modelInPhase1_.set(models.set(i, nullptr).ptr());
        }
        else
        {
            modelInPhase2_.set(models.set(i, nullptr).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ModelType>
Foam::SidedInterfacialModel<ModelType>::~SidedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ModelType>
const Foam::phaseInterface&
Foam::SidedInterfacialModel<ModelType>::interface() const
{
    return interface_;
}


template<class ModelType>
bool Foam::SidedInterfacialModel<ModelType>::haveModelInThe
(
    const phaseModel& phase
) const
{
    return
        interface_.index(phase) == 0
      ? modelInPhase1_.valid()
      : modelInPhase2_.valid();
}


template<class ModelType>
const ModelType& Foam::SidedInterfacialModel<ModelType>::modelInThe
(
    const phaseModel& phase
) const
{
    if (!haveModelInThe(phase))
    {
        FatalErrorInFunction
            << "There is no " << type() << " active for the "
            << phase.name() << " side of the "
            << interface_.name() << " interface"
            << exit(FatalError);
    }

    return
        interface_.index(phase) == 0
      ? modelInPhase1_()
      : modelInPhase2_();
}


template<class ModelType>
ModelType& Foam::SidedInterfacialModel<ModelType>::modelInThe
(
    const phaseModel& phase
)
{
    if (!haveModelInThe(phase))
    {
        FatalErrorInFunction
            << "There is no " << type() << " active for the "
            << phase.name() << " side of the "
            << interface_.name() << " interface"
            << exit(FatalError);
    }

    return
        interface_.index(phase) == 0
      ? modelInPhase1_()
      : modelInPhase2_();
}


template<class ModelType>
bool Foam::SidedInterfacialModel<ModelType>::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
