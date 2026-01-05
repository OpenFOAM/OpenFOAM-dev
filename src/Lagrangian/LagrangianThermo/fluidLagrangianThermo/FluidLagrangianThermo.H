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
    Foam::FluidLagrangianThermo

Description
    Fluid Lagrangian thermodynamic model implementation

SourceFiles
    FluidLagrangianThermo.C

\*---------------------------------------------------------------------------*/

#ifndef FluidLagrangianThermo_H
#define FluidLagrangianThermo_H

#include "BasicLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class FluidLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BaseThermo>
class FluidLagrangianThermo
:
    public BaseThermo
{
public:

    // Constructors

        //- Inherit constructors
        using BaseThermo::BaseThermo;


    //- Destructor
    virtual ~FluidLagrangianThermo();


    // Member Functions

        //- Update state
        virtual void correct(const LagrangianSubMesh&);


        // Thermodynamic State

            //- Compressibility for an injection [m^2/s^2]
            using BaseThermo::psi;


        // Derived Thermodynamic Properties

            //- Compressibility for an injection [m^2/s^2]
            virtual tmp<LagrangianSubScalarField> psi
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;


        // Transport State

            //- Dynamic viscosity [kg/m/s]
            using BaseThermo::mu;


        // Derived Transport Properties

            //- Dynamic viscosity for an injection [W/m/K]
            virtual tmp<LagrangianSubScalarField> mu
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FluidLagrangianThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
