/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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
    Foam::multicomponentThermo

Description
    Base-class for multi-component thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    multicomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentThermo_H
#define multicomponentThermo_H

#include "basicThermo.H"
#include "MulticomponentThermo.H"
#include "speciesTable.H"
#include "GeometricFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class multicomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentThermo
:
    virtual public basicThermo
{
public:

    //- Runtime type information
    TypeName("multicomponentThermo");


    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Destructor
    virtual ~multicomponentThermo();


    // Member Functions

        // Species set

            //- The table of species
            virtual const speciesTable& species() const = 0;

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const;

            //- The index of the default specie
            virtual label defaultSpecie() const = 0;


        // Active/Inactive handling

            //- Access the specie active flags
            virtual const List<bool>& speciesActive() const = 0;

            //- Synchronise the specie active flags
            void syncSpeciesActive() const;

            //- Set specie active
            inline void setSpecieActive(const label speciei) const;

            //- Set specie inactive
            inline void setSpecieInactive(const label speciei) const;

            //- Should the given specie be solved for? I.e., is it active and
            //  not the default specie?
            inline bool solveSpecie(const label speciei) const;


        // Mass fractions

            //- Access the mass-fraction fields
            virtual PtrList<volScalarField>& Y() = 0;

            //- Access the mass-fraction fields
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Access the mass-fraction field for a specie given by index
            inline volScalarField& Y(const label speciei);

            //- Access the mass-fraction field for a specie given by index
            inline const volScalarField& Y(const label speciei) const;

            //- Access the mass-fraction field for a specie given by name
            inline volScalarField& Y(const word& specieName);

            //- Access the mass-fraction field for a specie given by name
            inline const volScalarField& Y(const word& specieName) const;

            //- Access the specie index of the given mass-fraction field
            inline label specieIndex(const volScalarField& Yi) const;

            //- Normalise the mass fractions by clipping positive and deriving
            //  the default specie mass fraction from the other species.
            void normaliseY();


        // Specie molecular properties

            //- Molecular weight [kg/kmol]
            virtual scalar WiValue(const label speciei) const = 0;

            //- Molecular weight [kg/kmol]
            virtual dimensionedScalar Wi(const label speciei) const = 0;


        // Specie thermodynamic properties

            //- Density [kg/m^3]
            virtual scalar rhoi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rhoi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cpi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cpi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar hei
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> hei
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField::Internal> hei
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> hei
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar hsi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> hsi
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField::Internal> hsi
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hsi
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar hai
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<scalarField> hai
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField::Internal> hai
            (
                const label speciei,
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> hai
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual scalar hfiValue(const label speciei) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual dimensionedScalar hfi(const label speciei) const = 0;


        // Specie transport properties

            //- Thermal conductivity [W/m/K]
            virtual scalar kappai
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappai
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;
};


/*---------------------------------------------------------------------------*\
                     Class multicomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentThermo::implementation
:
    virtual public multicomponentThermo
{
protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- The name of the default specie. The mass fraction of which is
        //  derived from the other species rather than solved. Also used as the
        //  carrier specie in multicomponent diffusion.
        word defaultSpecieName_;

        //- The index of the default specie
        label defaultSpeciei_;

        //- List of specie active flags
        List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


    // Protected Member Functions

        //- Scale the mass fractions to sum to 1
        void correctMassFractions();


public:

    // Constructors

        //- Construct from dictionary, specie names, mesh and phase name
        implementation
        (
            const dictionary&,
            const wordList&,
            const fvMesh&,
            const word&
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- The table of species
        virtual const speciesTable& species() const;

        //- The index of the default specie
        virtual label defaultSpecie() const;

        //- Access the specie active flags
        virtual const List<bool>& speciesActive() const;

        //- Access the mass-fraction fields
        virtual PtrList<volScalarField>& Y();

        //- Access the mass-fraction fields
        virtual const PtrList<volScalarField>& Y() const;

        //- Get the slicer
        inline volScalarFieldListSlicer Yslicer() const;

        //- Get the composition of an internal cell
        inline scalarFieldListSlice cellComposition
        (
            const volScalarFieldListSlicer& Yslicer,
            const label celli
        ) const;

        //- Get the composition of a boundary face
        inline scalarFieldListSlice patchFaceComposition
        (
            const volScalarFieldListSlicer& Yslicer,
            const label patchi,
            const label facei
        ) const;

        //- Get the slicer for the given source
        inline PtrList<scalarField> Yslicer(const fvSource& source) const;

        //- Get the composition of a source cell
        inline scalarFieldListSlice sourceCellComposition
        (
            const PtrList<scalarField>& Yslicer,
            const label i
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
