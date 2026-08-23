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
    Foam::surfaceThermo

Description
    Base-class for Lagrangian surface thermodynamic models. This system holds a
    reference to a far-field thermo (either that of the carrier, or the thermo
    within the Lagrangian elements) and by default forwards properties from
    that. It also provides methods with which Lagrangian models can modify the
    surface thermodynamic state (temperature, pressure and species fractions).
    If the surface state is modified, the model then re-evaluates and caches
    properties at that state and returns those instead of the far-field values.

    This system also constructs a surface function by which the surface and
    far-field transport properties can be combined to produce transport
    properties which are representative of transfer processes (drag, heat
    transfer, phase change, ...) occurring between the surface and the far
    field. These are referred to as 'near' surface transport properties, and
    this class provides methods with which to access them also.

SourceFiles
    surfaceThermo.C
    surfaceThermoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceThermo_H
#define surfaceThermo_H

#include "basicLagrangianSubThermo.H"
#include "surfaceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianModel;

/*---------------------------------------------------------------------------*\
                      Class surfaceThermo Declaration
\*---------------------------------------------------------------------------*/

class surfaceThermo
:
    virtual public basicLagrangianSubThermo
{
protected:

    // Protected Member Data

        //- The sub-mesh index associated with any cached state
        mutable uint64_t subMeshIndex_;

        //- Has we deviated from the far thermo model?
        mutable bool surface_;


    // Protected Member Functions

        //- Reset the surface thermo if the sub-mesh has changed
        void setSubMesh(const LagrangianSubMesh& subMesh) const;

        //- Reset to the far thermo
        virtual void reset() const;

        //- Reset to the far thermo
        virtual void resetT() const = 0;

        //- Reset a derived property following a change in surface state
        template<template<class> class PrimitiveField>
        void resetDerived
        (
            tmp<LagrangianSubField<scalar, PrimitiveField>>& tf
        ) const;

        //- Reset derived properties following a change in surface state
        virtual void resetDerived() const;

        //- Reset derived properties following a change in surface state
        virtual void resetDerivedT() const = 0;

        //- Correct and cache derived properties for the given sub-mesh
        virtual void correct(const LagrangianSubMesh&) const = 0;

        //- Retrieve a derived property, correcting and caching if needed
        template<template<class> class PrimitiveField>
        const LagrangianSubField<scalar, PrimitiveField>& getDerived
        (
            const LagrangianSubMesh& subMesh,
            tmp<LagrangianSubField<scalar, PrimitiveField>>& tf
        ) const;


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    // Public Typedefs

        //- The derived name
        static word derivedThermoName()
        {
            return "surfaceThermo";
        }


    //- Runtime type information
    TypeName("surfaceThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceThermo,
        cloud,
        (const basicLagrangianSubThermo& farThermo),
        (farThermo)
    );


    // Constructors

        //- Default construct
        surfaceThermo();


    // Selectors

        //- Select a surface thermo
        static autoPtr<surfaceThermo> New
        (
            const basicLagrangianSubThermo& farThermo
        );

        //- Select a surface thermo
        template<class SurfaceThermo>
        static autoPtr<SurfaceThermo> New
        (
            const basicLagrangianSubThermo& farThermo
        );


    //- Destructor
    virtual ~surfaceThermo();


    // Member Functions

        //- Access the surface function
        virtual const Foam::surfaceFunction& surfaceFunction() const = 0;


        // Thermodynamic State

            //- Allow a model to set the temperature [K]
            virtual void setT
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarSubField>& tT
            ) = 0;

            //- Allow a model to set the temperature [K]
            virtual void setT
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarField>& tT
            );


        // Transport Properties

            //- Near thermal conductivity for a sub-mesh [W/m/K]
            virtual tmp<LagrangianSubScalarSubField> kappaNear
            (
                const LagrangianSubMesh& subMesh
            ) const = 0;
};

/*---------------------------------------------------------------------------*\
             Class surfaceThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class surfaceThermo::implementation
:
    virtual public surfaceThermo
{
protected:

    // Protected Data

        //- Bulk thermo
        const basicLagrangianSubThermo& farThermo_;

        //- Cached surface molecular weight [kg/kmol]
        mutable tmp<LagrangianSubScalarField> tWs_;

        //- Cached surface temperature [K]
        mutable tmp<LagrangianSubScalarSubField> tTs_;

        //- The model that set the surface temperature
        mutable const LagrangianModel* TmodelPtr_;

        //- Cached surface density [kg/m^3]
        mutable tmp<LagrangianSubScalarSubField> trhos_;

        //- Cached surface internal energy [J/kg]
        mutable tmp<LagrangianSubScalarSubField> tes_;

        //- Cached surface sensible enthalpy [J/kg]
        mutable tmp<LagrangianSubScalarField> thss_;

        //- Cached surface absolute enthalpy [J/kg]
        mutable tmp<LagrangianSubScalarField> thas_;

        //- Cached surface heat capacity at constant pressure [J/kg/K]
        mutable tmp<LagrangianSubScalarField> tCps_;

        //- Cached surface heat capacity at constant volume [J/kg/K]
        mutable tmp<LagrangianSubScalarSubField> tCvs_;

        //- Cached surface coefficient of thermal expansion [1/K]
        mutable tmp<LagrangianSubScalarField> talphavs_;

        //- Cached thermal conductivity [W/m/K]
        mutable tmp<LagrangianSubScalarSubField> tkappas_;

        //- Surface function
        mutable autoPtr<Foam::surfaceFunction> surfaceFunctionPtr_;


    // Protected Member Functions

        //- Reset to the far thermo
        virtual void resetT() const;

        //- Reset derived properties following a change in surface state
        virtual void resetDerivedT() const;


public:

    // Constructors

        //- Construct from a reference to the far thermo
        implementation(const basicLagrangianSubThermo& farThermo);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const;

        //- Access the mesh
        virtual const LagrangianMesh& mesh() const;

        //- Phase name
        virtual const word& phaseName() const;

        //- Access the surface function
        const Foam::surfaceFunction& surfaceFunction() const;


        // Molecular Properties

            //- Molecular weight for a sub-mesh [kg/kmol]
            virtual tmp<LagrangianSubScalarField> W
            (
                const LagrangianSubMesh&
            ) const;


        // Thermodynamic State

            //- Temperature for a sub-mesh [K]
            virtual tmp<LagrangianSubScalarSubField> T
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Allow a model to set the temperature [K]
            virtual void setT
            (
                const LagrangianModel& model,
                const tmp<LagrangianSubScalarSubField>& tT
            );


        // Thermodynamic Properties

            //- Density for a sub-mesh [kg/m^3]
            virtual tmp<LagrangianSubScalarSubField> rho
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Internal energy for a sub-mesh [J/kg]
            virtual tmp<LagrangianSubScalarSubField> e
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Sensible enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> hs
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Absolute enthalpy for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> ha
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Heat capacity at constant pressure for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cp
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Heat capacity at constant volume for a sub-mesh [J/kg/K]
            virtual tmp<LagrangianSubScalarSubField> Cv
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Coefficient of thermal expansion for a sub-mesh [1/K]
            virtual tmp<LagrangianSubScalarField> alphav
            (
                const LagrangianSubMesh& subMesh
            ) const;


        // Transport Properties

            //- Thermal conductivity for a sub-mesh [W/m/K]
            virtual tmp<LagrangianSubScalarSubField> kappa
            (
                const LagrangianSubMesh& subMesh
            ) const;

            //- Near thermal conductivity for a sub-mesh [W/m/K]
            virtual tmp<LagrangianSubScalarSubField> kappaNear
            (
                const LagrangianSubMesh& subMesh
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
