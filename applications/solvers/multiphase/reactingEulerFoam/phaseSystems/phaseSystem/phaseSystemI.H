/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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


inline const Foam::phaseSystem::phasePairTable&
Foam::phaseSystem::phasePairs() const
{
    return phasePairs_;
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


inline Foam::fv::IOoptionList& Foam::phaseSystem::fvOptions() const
{
    return fvOptions_;
}


// ************************************************************************* //
