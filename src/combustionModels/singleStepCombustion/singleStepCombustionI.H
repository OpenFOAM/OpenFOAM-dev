/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

template<class ThermoType>
const multiComponentMixture<ThermoType>&
singleStepCombustion<ThermoType>::mixture() const
{
    return mixture_;
}


template<class ThermoType>
inline const Foam::Reaction<ThermoType>&
singleStepCombustion<ThermoType>::reaction() const
{
    return reaction_();
}


template<class ThermoType>
inline const dimensionedScalar&
singleStepCombustion<ThermoType>::stoicRatio() const
{
    return stoicRatio_;
}


template<class ThermoType>
inline const dimensionedScalar& singleStepCombustion<ThermoType>::s() const
{
    return s_;
}


template<class ThermoType>
inline const dimensionedScalar& singleStepCombustion<ThermoType>::qFuel() const
{
    return qFuel_;
}


template<class ThermoType>
inline const List<scalar>&
singleStepCombustion<ThermoType>::specieStoichCoeffs() const
{
    return specieStoichCoeffs_;
}


template<class ThermoType>
inline tmp<volScalarField> singleStepCombustion<ThermoType>::fres
(
    const label index
) const
{
    return fres_[index];
}


template<class ThermoType>
inline label singleStepCombustion<ThermoType>::inertIndex() const
{
    return inertIndex_;
}


template<class ThermoType>
inline label singleStepCombustion<ThermoType>::fuelIndex() const
{
    return fuelIndex_;
}


template<class ThermoType>
inline const List<int>& singleStepCombustion<ThermoType>::specieProd() const
{
    return specieProd_;
}


template<class ThermoType>
inline const scalarList& singleStepCombustion<ThermoType>::Yprod0() const
{
    return Yprod0_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// ************************************************************************* //
