/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::psiuMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    psiuMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiuMulticomponentThermo_H
#define psiuMulticomponentThermo_H

#include "psiThermo.H"
#include "PsiuMulticomponentThermo.H"
#include "speciesTable.H"
#include "GeometricFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiuMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo
:
    virtual public psiThermo
{
protected:

    // Protected Member Functions

        //- Return the unburnt enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList heuBoundaryTypes();

        //- ...
        void heuBoundaryCorrection(volScalarField& heu);


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
            PsiuMulticomponentThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName;


    //- Runtime type information
    TypeName("psiuMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            psiuMulticomponentThermo,
            fvMesh,
            (const fvMesh& mesh, const word& phaseName),
            (mesh, phaseName)
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiuMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiuMulticomponentThermo();


    // Member Functions

        // Species set

            //- The table of species
            virtual const speciesTable& species() const = 0;

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const;


        // Mass fractions

            //- Access the mass-fraction fields
            virtual PtrList<volScalarField>& Y() = 0;

            //- Access the mass-fraction fields
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return the mass-fraction field for a specie given by name
            inline volScalarField& Y(const word& specieName);

            //- Return the const mass-fraction field for a specie given by name
            inline const volScalarField& Y(const word& specieName) const;

            //- Return the residual fraction of fuel in the burnt mixture
            virtual tmp<volScalarField> fres() const = 0;

            //- Reset the mixture to an unburnt state and update EGR
            virtual void reset() = 0;


        // Thermodynamic state

            //- Unburnt gas temperature [K]
            virtual const volScalarField& Tu() const = 0;

            //- Unburnt gas enthalpy [J/kg]
            virtual const volScalarField& heu() const = 0;

            //- Unburnt gas enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heu() = 0;


        // Derived thermodynamic properties

            //- Unburnt gas enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Unburnt gas enthalpy for patch [J/kg]
            virtual tmp<scalarField> heu
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Standard enthalpy of reaction [J/kg]
            virtual tmp<volScalarField> hr() const = 0;

            //- Burnt gas temperature [K]
            virtual tmp<volScalarField> Tb() const = 0;

            //- Unburnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psiu() const = 0;

            //- Burnt gas compressibility [s^2/m^2]
            virtual tmp<volScalarField> psib() const = 0;

            //- Unburnt gas density [kg/m^3]
            tmp<volScalarField> rhou() const;

            //- Burnt gas density [kg/m^3]
            tmp<volScalarField> rhob() const;


        // Derived transport properties

            //- Dynamic viscosity of unburnt gas [kg/m/s]
            virtual tmp<volScalarField> muu() const = 0;

            //- Dynamic viscosity of burnt gas [kg/m/s]
            virtual tmp<volScalarField> mub() const = 0;
};


/*---------------------------------------------------------------------------*\
           Class psiuMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo::implementation
:
    virtual public psiuMulticomponentThermo
{
protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


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


/*---------------------------------------------------------------------------*\
                Class psiuMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiuMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiThermo::implementation,
    public psiuMulticomponentThermo::implementation
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName),
            psiThermo::implementation(dict, mesh, phaseName),
            psiuMulticomponentThermo::implementation
            (
                dict,
                MixtureType::specieNames(),
                mesh,
                phaseName
            )
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "psiuMulticomponentThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
