/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "singleStepReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::reaction& singleStepReaction::singleReaction() const
{
    return reaction_;
}


inline const dimensionedScalar& singleStepReaction::stoicRatio() const
{
    return stoicRatio_;
}


inline const dimensionedScalar& singleStepReaction::s() const
{
    return s_;
}


inline const dimensionedScalar& singleStepReaction::qFuel() const
{
    return qFuel_;
}


inline const List<scalar>& singleStepReaction::specieStoichCoeffs() const
{
    return specieStoichCoeffs_;
}


inline tmp<volScalarField> singleStepReaction::fres
(
    const label index
) const
{
    return fres_[index];
}


inline label singleStepReaction::fuelIndex() const
{
    return fuelIndex_;
}


inline const List<int>& singleStepReaction::specieProd() const
{
    return specieProd_;
}


inline const scalarList& singleStepReaction::Yprod0() const
{
    return Yprod0_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace Foam

// ************************************************************************* //
