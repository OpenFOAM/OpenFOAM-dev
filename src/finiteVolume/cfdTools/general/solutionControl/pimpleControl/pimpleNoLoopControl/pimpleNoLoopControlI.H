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

#include "pimpleLoop.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::pimpleNoLoopControl::moveMeshOuterCorrectors() const
{
    return moveMeshOuterCorrectors_;
}


inline bool Foam::pimpleNoLoopControl::simpleRho() const
{
    return simpleRho_;
}


inline bool Foam::pimpleNoLoopControl::predictTransport() const
{
    return !transportPredictionFirst_ || firstIter();
}


inline bool Foam::pimpleNoLoopControl::correctTransport() const
{
    return !transportCorrectionFinal_ || finalIter();
}


inline bool Foam::pimpleNoLoopControl::firstIter() const
{
    return pimpleLoopPtr_ ? pimpleLoopPtr_->firstIter() : true;
}


inline bool Foam::pimpleNoLoopControl::finalIter() const
{
    return pimpleLoopPtr_ ? pimpleLoopPtr_->finalIter() : true;
}


// ************************************************************************* //
