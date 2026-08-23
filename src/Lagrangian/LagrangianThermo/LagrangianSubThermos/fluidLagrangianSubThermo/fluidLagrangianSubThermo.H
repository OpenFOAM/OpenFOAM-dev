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
    Foam::fluidLagrangianThermo

Description
    Abstract interface which provides access to fluid thermodynamic properties
    on Lagrangian sub-meshes.

SourceFiles
    fluidLagrangianSubThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidLagrangianSubThermo_H
#define fluidLagrangianSubThermo_H

#include "basicLagrangianSubThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fluidLagrangianSubThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidLagrangianSubThermo
:
    virtual public basicLagrangianSubThermo
{
public:

    //- Runtime type information
    TypeName("fluidLagrangianSubThermo");


    //- Destructor
    virtual ~fluidLagrangianSubThermo();


    // Member Functions

        // Thermodynamic State

            //- Pressure for a sub-mesh [Pa]
            virtual tmp<LagrangianSubScalarSubField> p
            (
                const LagrangianSubMesh&
            ) const = 0;


        // Thermodynamic Properties

            //- Compressibility for a sub-mesh [m^2/s^2]
            virtual tmp<LagrangianSubScalarSubField> psi
            (
                const LagrangianSubMesh&
            ) const = 0;


        // Transport Properties

            //- Dynamic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarSubField> mu
            (
                const LagrangianSubMesh&
            ) const = 0;

            //- Kinematic viscosity for a sub-mesh [m^2/s]
            virtual tmp<LagrangianSubScalarField> nu
            (
                const LagrangianSubMesh&
            ) const;

            //- Prandtl number for a sub-mesh []
            virtual tmp<LagrangianSubScalarField> Pr
            (
                const LagrangianSubMesh&
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
