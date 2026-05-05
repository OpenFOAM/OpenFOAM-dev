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
    Foam::multicomponentLagrangianThermo

Description
    Base-class for multicomponent Lagrangian thermodynamic models

SourceFiles
    multicomponentLagrangianThermoI.H
    multicomponentLagrangianThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentLagrangianThermo_H
#define multicomponentLagrangianThermo_H

#include "basicLagrangianThermo.H"
#include "MulticomponentLagrangianThermo.H"
#include "speciesTable.H"
#include "DimensionedFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multicomponentLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentLagrangianThermo
:
    virtual public basicLagrangianThermo
{
public:

    //- Runtime type information
    TypeName("multicomponentLagrangianThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        multicomponentLagrangianThermo,
        LagrangianMesh,
        (const LagrangianMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    // Selectors

        //- Select a multicomponent thermo
        static autoPtr<multicomponentLagrangianThermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~multicomponentLagrangianThermo();


    // Member Functions

        // Species set

            //- Return the table of species
            virtual const speciesTable& species() const = 0;

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const;

            //- The index of the default specie
            virtual label defaultSpecie() const = 0;


        // Mass fractions

            //- Access the mass-fraction fields
            virtual PtrList<LagrangianScalarDynamicField>& Y() = 0;

            //- Access the mass-fraction fields
            virtual const PtrList<LagrangianScalarDynamicField>& Y() const = 0;

            //- Access the mass-fraction field for a specie given by index
            inline LagrangianScalarDynamicField& Y(const label);

            //- Access the mass-fraction field for a specie given by index
            inline const LagrangianScalarDynamicField& Y(const label) const;

            //- Access the mass-fraction field for a specie given by name
            inline LagrangianScalarDynamicField& Y(const word&);

            //- Access the mass-fraction field for a specie given by name
            inline const LagrangianScalarDynamicField& Y(const word&) const;

            //- Access the specie index of the given mass-fraction field
            inline label specieIndex(const LagrangianScalarDynamicField&) const;

            //- Normalise the mass fractions within the given sub-mesh by
            //  clipping positive and deriving the default specie mass fraction
            //  from the other species
            virtual void normaliseY(const LagrangianSubMesh&) = 0;


        // Specie molecular properties

            //- Molecular weight [kg/kmol]
            virtual dimensionedScalar Wi(const label speciei) const = 0;


        // Specie thermodynamic properties

            //- Density [kg/m^3]
            virtual tmp<LagrangianSubScalarField> rhoi
            (
                const label speciei,
                const LagrangianSubScalarField& p,
                const LagrangianSubScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<LagrangianSubScalarField> hsi
            (
                const label speciei,
                const LagrangianSubScalarField& p,
                const LagrangianSubScalarField& T
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual dimensionedScalar hfi(const label speciei) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cvi
            (
                const label speciei,
                const LagrangianSubScalarField& p,
                const LagrangianSubScalarField& T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cpi
            (
                const label speciei,
                const LagrangianSubScalarField& p,
                const LagrangianSubScalarField& T
            ) const = 0;
};


/*---------------------------------------------------------------------------*\
        Class multicomponentLagrangianThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class multicomponentLagrangianThermo::implementation
:
    virtual public multicomponentLagrangianThermo
{
protected:

    // Protected data

        //- The name of the default specie. The mass fraction of which is
        //  derived from the other species rather than solved. Also used as the
        //  carrier specie in multicomponent diffusion.
        word defaultSpecieName_;

        //- The index of the default specie
        label defaultSpeciei_;

        //- Species mass fractions
        PtrList<LagrangianScalarDynamicField> Y_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation
        (
            const dictionary&,
            const speciesTable&,
            const LagrangianMesh&,
            const word&
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Species set

            //- The index of the default specie
            virtual label defaultSpecie() const;


        // Mass fractions

            //- Access the mass-fraction fields
            virtual PtrList<LagrangianScalarDynamicField>& Y();

            //- Access the mass-fraction fields
            virtual const PtrList<LagrangianScalarDynamicField>& Y() const;

            //- Normalise the mass fractions within the given sub-mesh by
            //  clipping positive and deriving the default specie mass fraction
            //  from the other species
            virtual void normaliseY(const LagrangianSubMesh&);


        // Composition

            //- Get the slicer
            inline DimensionedFieldListSlicer<scalar> Yslicer() const;

            //- Get the composition of an element
            inline scalarFieldListSlice elementComposition
            (
                const DimensionedFieldListSlicer<scalar>& Yslicer,
                const label i
            ) const;

            //- Get the slicer for an injection
            inline
            DimensionedFieldListAndSlicer<scalar, LagrangianSubMesh> Yslicer
            (
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const;

            //- Get the composition of an element in an injection
            inline scalarFieldListSlice injectionElementComposition
            (
                const DimensionedFieldListAndSlicer<scalar, LagrangianSubMesh>&,
                const label subi
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentLagrangianThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
