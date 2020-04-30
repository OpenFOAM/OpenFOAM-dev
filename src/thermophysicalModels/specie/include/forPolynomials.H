/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef forPolynomials_H
#define forPolynomials_H

#include "specie.H"

#include "icoPolynomial.H"

#include "eConstThermo.H"
#include "hConstThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "logPolynomialTransport.H"
#include "polynomialTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forPolynomialEqns(Mu, He, Cp, Macro, Args...)                          \
    forThermo(Mu, He, Cp, icoPolynomial, specie, Macro, Args);

#define forPolynomialEnergiesAndThermos(Mu, Macro, Args...)                    \
    forPolynomialEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);        \
    forPolynomialEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);  \
    forPolynomialEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args)

#define forPolynomialTransports(Macro, Args...)                                \
    forPolynomialEnergiesAndThermos(logPolynomialTransport, Macro, Args);      \
    forPolynomialEnergiesAndThermos(polynomialTransport, Macro, Args)

#define forPolynomials(Macro, Args...)                                         \
    forPolynomialTransports(Macro, Args)

#endif

// ************************************************************************* //
