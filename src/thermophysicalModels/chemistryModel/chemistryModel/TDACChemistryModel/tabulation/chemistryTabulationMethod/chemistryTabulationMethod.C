/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "chemistryTabulationMethod.H"
#include "TDACChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryTabulationMethod<CompType, ThermoType>::chemistryTabulationMethod
(
    const dictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    dict_(dict),
    coeffsDict_(dict.subDict("tabulation")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    log_(coeffsDict_.lookupOrDefault<Switch>("log", false)),
    chemistry_(chemistry),
    tolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance", 1e-4))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::chemistryTabulationMethod<CompType, ThermoType>::
~chemistryTabulationMethod()
{}


// ************************************************************************* //
