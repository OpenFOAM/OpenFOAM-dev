/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

inline const Foam::fvMesh& Foam::twoPhaseSystem::mesh() const
{
    return mesh_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::phase1() const
{
    return phase1_;
}


inline Foam::phaseModel& Foam::twoPhaseSystem::phase1()
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::phase2() const
{
    return phase2_;
}


inline Foam::phaseModel& Foam::twoPhaseSystem::phase2()
{
    return phase2_;
}


inline const Foam::phaseModel& Foam::twoPhaseSystem::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase == &phase1_)
    {
        return phase2_;
    }
    else
    {
        return phase1_;
    }
}


inline const Foam::surfaceScalarField& Foam::twoPhaseSystem::phi() const
{
    return phi_;
}


inline Foam::surfaceScalarField& Foam::twoPhaseSystem::phi()
{
    return phi_;
}


inline const Foam::volScalarField& Foam::twoPhaseSystem::dgdt() const
{
    return dgdt_;
}


inline Foam::volScalarField& Foam::twoPhaseSystem::dgdt()
{
    return dgdt_;
}


inline Foam::tmp<Foam::surfaceScalarField>& Foam::twoPhaseSystem::pPrimeByA()
{
    return pPrimeByA_;
}


// ************************************************************************* //
