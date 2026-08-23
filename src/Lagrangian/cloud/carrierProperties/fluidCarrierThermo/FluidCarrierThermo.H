/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::FluidCarrierThermo

Description
    Thermodynamic model of the carrier phase of a Lagrangian cloud.

SourceFiles
    FluidCarrierThermo.C

\*---------------------------------------------------------------------------*/

#ifndef FluidCarrierThermo_H
#define FluidCarrierThermo_H

#include "carried.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FluidCarrierThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class FluidCarrierThermo
:
    public BaseThermo
{
public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        FluidCarrierThermo
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const basicThermo& thermo
        );


    //- Destructor
    virtual ~FluidCarrierThermo();


    // Member Functions

        // Thermodynamic Properties

            //- Compressibility for a sub-mesh [m^2/s^2]
            virtual tmp<LagrangianSubScalarSubField> psi
            (
                const LagrangianSubMesh&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FluidCarrierThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
