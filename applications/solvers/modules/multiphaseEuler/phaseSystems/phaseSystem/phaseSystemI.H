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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}


inline const Foam::pimpleNoLoopControl& Foam::phaseSystem::pimple() const
{
    return pimple_;
}


inline const Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases() const
{
    return phaseModels_;
}


inline Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases()
{
    return phaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases() const
{
    return movingPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases()
{
    return movingPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases() const
{
    return stationaryPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases()
{
    return stationaryPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases() const
{
    return anisothermalPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases()
{
    return anisothermalPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multicomponentPhases() const
{
    return multicomponentPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multicomponentPhases()
{
    return multicomponentPhaseModels_;
}


inline const Foam::phaseModel& Foam::phaseSystem::otherPhase
(
    const phaseModel& phase
) const
{
    if (phaseModels_.size() != 2)
    {
        FatalErrorInFunction
            << "Call from a two-phase model in a multi-phase system."
            << exit(FatalError);
    }

    if (&phase == &phaseModels_[0])
    {
        return phaseModels_[1];
    }
    else
    {
        return phaseModels_[0];
    }
}


inline const Foam::surfaceScalarField& Foam::phaseSystem::phi() const
{
    return phi_;
}


inline Foam::surfaceScalarField& Foam::phaseSystem::phi()
{
    return phi_;
}


inline const Foam::volScalarField& Foam::phaseSystem::dpdt() const
{
    return dpdt_;
}


inline Foam::volScalarField& Foam::phaseSystem::dpdt()
{
    return dpdt_;
}


inline const Foam::IOMRFZoneList& Foam::phaseSystem::MRF() const
{
    return MRF_;
}


inline Foam::fvModels& Foam::phaseSystem::fvModels(fvMesh& mesh)
{
    return Foam::fvModels::New(mesh);
}


inline const Foam::fvModels& Foam::phaseSystem::fvModels() const
{
    return Foam::fvModels::New(mesh_);
}


inline Foam::fvConstraints& Foam::phaseSystem::fvConstraints(fvMesh& mesh)
{
    return Foam::fvConstraints::New(mesh);
}


inline const Foam::fvConstraints& Foam::phaseSystem::fvConstraints() const
{
    return Foam::fvConstraints::New(mesh_);
}


inline const Foam::dimensionedScalar& Foam::phaseSystem::deltaN() const
{
    return deltaN_;
}


// ************************************************************************* //
