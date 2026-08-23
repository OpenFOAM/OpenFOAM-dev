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
    Foam::surfaceFunction

Description
    Functions which combine values of transport properties evaluated on the
    surface and in the far-field to generate a representative value appropriate
    for calculating transfer processes.

SourceFiles
    surfaceFunction.C
    surfaceFunctionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFunction_H
#define surfaceFunction_H

#include "runTimeSelectionTables.H"
#include "LagrangianMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFunction Declaration
\*---------------------------------------------------------------------------*/

class surfaceFunction
{
    // Private Static Member Functions

        //- Convert a temporary sub-field in to a given primitive field type.
        //  Either forwards if PrimitiveField is Field, or uses toSubField if
        //  PrimitiveField is SubField.
        template<template<class> class PrimitiveField>
        static inline tmp<LagrangianSubField<scalar, PrimitiveField>>
        toPrimitiveField
        (
            const tmp<LagrangianSubScalarField>&
        );


    // Private Member Functions

        //- Return coefficients from which to construct a near-surface thermal
        //  conductivity
        virtual void thermalConductivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const = 0;

        //- Return coefficients from which to construct a near-surface dynamic
        //  diffusivity or viscosity
        virtual void dynamicDiffusivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const = 0;

        //- Return coefficients from which to construct a near-surface kinematic
        //  diffusivity or viscosity
        virtual void kinematicDiffusivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const = 0;

        //- Return a near-surface transport property, given the far and
        //  surface transport properties
        template<class Thermo, template<class> class PrimitiveField>
        tmp<LagrangianSubField<scalar, PrimitiveField>> transportProperty
        (
            const LagrangianSubMesh& subMesh,
            void (surfaceFunction::*coeffs)
            (
                const LagrangianSubMesh& subMesh,
                autoPtr<dimensionedScalar>&,
                tmp<LagrangianSubScalarField>&
            ) const,
            const Thermo& farThermo,
            const Thermo& surfaceThermo,
            tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
            (
                const LagrangianSubMesh& subMesh
            ) const
        ) const;

        //- Return a near-surface transport property, given the far and
        //  surface transport properties
        template<class FarFunction, class SurfaceFunction>
        tmp<LagrangianSubScalarField> transportProperty
        (
            const LagrangianSubMesh& subMesh,
            void (surfaceFunction::*coeffs)
            (
                const LagrangianSubMesh& subMesh,
                autoPtr<dimensionedScalar>&,
                tmp<LagrangianSubScalarField>&
            ) const,
            const FarFunction& farFunction,
            const SurfaceFunction& surfaceFunction
        ) const;


protected:

    // Protected Member Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


public:

    //- Runtime type information
    TypeName("surfaceFunction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceFunction,
            word,
            (const LagrangianMesh& mesh),
            (mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceFunction,
            dictionary,
            (const LagrangianMesh& mesh, const dictionary& dict),
            (mesh, dict)
        );


    // Constructors

        //- Construct for a given mesh
        surfaceFunction(const LagrangianMesh& mesh);

        //- Disallow default bitwise copy construction
        surfaceFunction(const surfaceFunction&) = delete;


    // Selectors

        //- Select given a mesh and a dictionary
        static autoPtr<surfaceFunction> New
        (
            const LagrangianMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceFunction();


    // Member Functions

        //- Return a near-surface thermal conductivity, given the far and
        //  surface thermal conductivities
        template<class Thermo, template<class> class PrimitiveField>
        tmp<LagrangianSubField<scalar, PrimitiveField>> thermalConductivity
        (
            const LagrangianSubMesh& subMesh,
            const Thermo& farThermo,
            const Thermo& surfaceThermo,
            tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
            (
                const LagrangianSubMesh& subMesh
            ) const
        ) const;

        //- Return a near-surface dynamic diffusivity, given the far and
        //  surface dynamic diffusivities
        template<class Thermo, template<class> class PrimitiveField>
        tmp<LagrangianSubField<scalar, PrimitiveField>> dynamicDiffusivity
        (
            const LagrangianSubMesh& subMesh,
            const Thermo& farThermo,
            const Thermo& surfaceThermo,
            tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
            (
                const LagrangianSubMesh& subMesh
            ) const
        ) const;

        //- Return a near-surface kinematic diffusivity, given the far and
        //  surface kinematic diffusivities
        template<class Thermo, template<class> class PrimitiveField>
        tmp<LagrangianSubField<scalar, PrimitiveField>> kinematicDiffusivity
        (
            const LagrangianSubMesh& subMesh,
            const Thermo& farThermo,
            const Thermo& surfaceThermo,
            tmp<LagrangianSubField<scalar, PrimitiveField>> (Thermo::*property)
            (
                const LagrangianSubMesh& subMesh
            ) const
        ) const;

        //- Return a near-surface kinematic diffusivity, given the far and
        //  surface kinematic diffusivities
        template<class FarFunction, class SurfaceFunction>
        tmp<LagrangianSubScalarField> kinematicDiffusivity
        (
            const LagrangianSubMesh& subMesh,
            const FarFunction& farFunction,
            const SurfaceFunction& surfaceFunction
        ) const;
};


//- Specialise forward of LagrangianSubScalarField
template<>
inline tmp<LagrangianSubScalarField>
surfaceFunction::toPrimitiveField<Field>
(
    const tmp<LagrangianSubScalarField>&
);

//- Specialise conversion from LagrangianSubScalarField to
//  LagrangianSubScalarSubFIeld
template<>
inline tmp<LagrangianSubScalarSubField>
surfaceFunction::toPrimitiveField<SubField>
(
    const tmp<LagrangianSubScalarField>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFunctionI.H"

#ifdef NoRepository
    #include "surfaceFunctionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
