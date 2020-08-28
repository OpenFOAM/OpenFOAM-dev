/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

Class
    Foam::sensibleEnthalpy

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleEnthalpy_H
#define sensibleEnthalpy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sensibleEnthalpy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensibleEnthalpy
{

public:

    // Constructors

        //- Construct
        sensibleEnthalpy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleEnthalpy";
        }

        // Fundamental properties

            static word energyName()
            {
                return "h";
            }

            // Heat capacity at constant pressure [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cp(p, T);
            }

            //- Cp/Cp []
            scalar CpByCpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return 1;
            }

            // Sensible enthalpy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Hs(p, T);
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar THE
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.THs(h, p, T0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
