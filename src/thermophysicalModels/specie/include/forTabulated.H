/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2025 OpenFOAM Foundation
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

#ifndef forTabulated_H
#define forTabulated_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specie.H"

#include "rhoTabulated.H"
#include "icoTabulated.H"

#include "hTabulatedThermo.H"
#include "eTabulatedThermo.H"
#include "hIcoTabulatedThermo.H"
#include "eIcoTabulatedThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "tabulatedTransport.H"
#include "icoTabulatedTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forTabulatedEqns(Mu, He, Cp, Macro, Args...)                           \
    forThermo(Mu, He, Cp, rhoTabulated, specie, Macro, Args)

#define forIcoTabulatedEqns(Mu, He, Cp, Macro, Args...)                        \
    forThermo(Mu, He, Cp, icoTabulated, specie, Macro, Args)

#define forTabulatedEnergiesAndThermos(Mu, Macro, Args...)                     \
    forTabulatedEqns(Mu, sensibleEnthalpy, hTabulatedThermo, Macro, Args);     \
    forTabulatedEqns(Mu, sensibleInternalEnergy, eTabulatedThermo, Macro, Args)

#define forIcoTabulatedEnergiesAndThermos(Mu, Args...)                         \
    forIcoTabulatedEqns(Mu, sensibleEnthalpy, hIcoTabulatedThermo, Args);      \
    forIcoTabulatedEqns(Mu, sensibleInternalEnergy, eIcoTabulatedThermo, Args)

#define forTabulatedTransport(Macro, Args...)                                  \
    forTabulatedEnergiesAndThermos(tabulatedTransport, Macro, Args)

#define forIcoTabulatedTransport(Macro, Args...)                               \
    forIcoTabulatedEnergiesAndThermos(icoTabulatedTransport, Macro, Args)

#define forTabulated(Macro, Args...)                                           \
    forTabulatedTransport(Macro, Args);                                        \
    forIcoTabulatedTransport(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
