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
    Foam::fluidSurfaceThermo

Description
    Base-class for fluid Lagrangian surface thermodynamic models

SourceFiles
    fluidSurfaceThermo.C
    fluidSurfaceThermos.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSurfaceThermo_H
#define fluidSurfaceThermo_H

#include "surfaceThermo.H"
#include "fluidLagrangianSubThermo.H"
#include "pureSurfaceThermo.H"
#include "SurfaceThermo.H"
#include "FluidSurfaceThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fluidSurfaceThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidSurfaceThermo
:
    virtual public surfaceThermo,
    virtual public fluidLagrangianSubThermo
{
protected:

    // Protected Member Functions

        //- Reset to the far thermo
        virtual void reset() const;

        //- Reset to the far thermo
        virtual void resetp() const = 0;

        //- Reset derived properties following a change in surface state
        virtual void resetDerived() const;

        //- Reset derived properties following a change in surface state
        virtual void resetDerivedp() const = 0;


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            FluidSurfaceThermo
            <
                SurfaceThermo
                <
                    MixtureType,
                    composite
                >
            >;


    //- Runtime type information
    TypeName("fluidSurfaceThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidSurfaceThermo,
        cloud,
        (const basicLagrangianSubThermo& farThermo),
        (farThermo)
    );


    // Selectors

        //- Select a surface thermo
        static autoPtr<fluidSurfaceThermo> New
        (
            const fluidLagrangianSubThermo& farFluidThermo
        );


    //- Destructor
    virtual ~fluidSurfaceThermo();


    // Member Functions

        // Thermodynamic State

            //- Allow a model to set the pressure [K]
            virtual void setp
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarSubField>& tp
            ) = 0;

            //- Allow a model to set the pressure [K]
            virtual void setp
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarField>& tp
            );


        // Transport Properties

            //- Near dynamic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarSubField> muNear
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Near kinematic viscosity for a sub-mesh [m^2/s]
            virtual tmp<LagrangianSubScalarField> nuNear
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;

            //- Near Prandtl number for a sub-mesh []
            virtual tmp<LagrangianSubScalarField> PrNear
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;
};

/*---------------------------------------------------------------------------*\
             Class fluidSurfaceThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidSurfaceThermo::implementation
:
    virtual public fluidSurfaceThermo
{
protected:

    // Protected Data

        //- Bulk thermo
        const fluidLagrangianSubThermo& farFluidThermo_;

        //- Surface pressure
        mutable tmp<LagrangianSubScalarSubField> tps_;

        //- The model that set the surface pressure
        mutable const LagrangianModel* pModelPtr_;

        //- Cached surface compressibility [m^2/s^2]
        mutable tmp<LagrangianSubScalarSubField> tpsis_;

        //- Cached surface dynamic viscosity [kg/m/s]
        mutable tmp<LagrangianSubScalarSubField> tmus_;


    // Protected Member Functions

        //- Reset to the far thermo
        virtual void resetp() const;

        //- Reset derived properties following a change in surface state
        virtual void resetDerivedp() const;

public:

    // Constructors

        //- Construct from a reference to the far thermo
        implementation(const basicLagrangianSubThermo& farThermo);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Thermodynamic State

            //- Pressure for a sub-mesh [Pa]
            virtual tmp<LagrangianSubScalarSubField> p
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Allow a model to set the temperature [K]
            virtual void setp
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarSubField>& tp
            );


        // Thermodynamic Properties

            //- Compressibility for a sub-mesh [m^2/s^2]
            virtual tmp<LagrangianSubScalarSubField> psi
            (
                const LagrangianSubMesh&
            ) const;


        // Transport Properties

            //- Dynamic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarSubField> mu
            (
                const LagrangianSubMesh&
            ) const;

            //- Near dynamic viscosity for a sub-mesh [kg/m/s]
            virtual tmp<LagrangianSubScalarSubField> muNear
            (
                const LagrangianSubMesh&
            ) const;

            //- Near kinematic viscosity for a sub-mesh [m^2/s]
            virtual tmp<LagrangianSubScalarField> nuNear
            (
                const LagrangianSubMesh&
            ) const;

            //- Near Prandtl number for a sub-mesh []
            virtual tmp<LagrangianSubScalarField> PrNear
            (
                const LagrangianSubMesh&
            ) const;
};


/*---------------------------------------------------------------------------*\
                 Class fluidSurfaceThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidSurfaceThermo::composite
:
    public surfaceThermo::implementation,
    public fluidSurfaceThermo::implementation,
    public pureSurfaceThermo
{
public:

    // Constructors

        //- Construct from a reference to the far thermo
        composite(const basicLagrangianSubThermo& farThermo)
        :
            surfaceThermo::implementation(farThermo),
            fluidSurfaceThermo::implementation(farThermo)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
