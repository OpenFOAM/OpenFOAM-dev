/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::reaction& singleStepCombustion::singleReaction() const
{
    return reaction_;
}


inline const dimensionedScalar& singleStepCombustion::stoicRatio() const
{
    return stoicRatio_;
}


inline const dimensionedScalar& singleStepCombustion::s() const
{
    return s_;
}


inline const dimensionedScalar& singleStepCombustion::qFuel() const
{
    return qFuel_;
}


inline const List<scalar>& singleStepCombustion::specieStoichCoeffs() const
{
    return specieStoichCoeffs_;
}


inline tmp<volScalarField> singleStepCombustion::fres
(
    const label index
) const
{
    return fres_[index];
}


inline label singleStepCombustion::fuelIndex() const
{
    return fuelIndex_;
}


inline const List<int>& singleStepCombustion::specieProd() const
{
    return specieProd_;
}


inline const scalarList& singleStepCombustion::Yprod0() const
{
    return Yprod0_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// ************************************************************************* //
