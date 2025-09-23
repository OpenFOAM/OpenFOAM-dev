/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::ubPsiMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    ubPsiMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ubPsiMulticomponentThermo_H
#define ubPsiMulticomponentThermo_H

#include "psiThermo.H"
#include "UBPsiMulticomponentThermo.H"
#include "speciesTable.H"
#include "DimensionedFieldListSlicer.H"
#include "GeometricFieldListSlicer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ubPsiMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class ubPsiMulticomponentThermo
:
    virtual public psiThermo
{

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
            UBPsiMulticomponentThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName()
        {
            return "ubPsiThermo";
        }


    //- Runtime type information
    TypeName("ubPsiMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            ubPsiMulticomponentThermo,
            fvMesh,
            (const fvMesh& mesh, const word& phaseName),
            (mesh, phaseName)
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<ubPsiMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~ubPsiMulticomponentThermo();


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

            //- Return the fuel-oxidant equivalence ratio
            virtual tmp<volScalarField> Phi() const = 0;

            //- Return the mixture-fraction
            virtual tmp<volScalarField> ft() const = 0;

            //- Return the fuel mass-fraction
            virtual tmp<volScalarField> fu() const = 0;

            //- Return the exhaust gas recirculation mixture-fraction
            virtual tmp<volScalarField> egr() const = 0;

            //- Reset the mixture to an unburnt state and update EGR
            virtual void reset
            (
                volScalarField& b,
                volScalarField& c,
                const PtrList<volScalarField>& Yb,
                const volScalarField& heb
            ) = 0;

        //- Standard enthalpy of formation [J/kg]
        virtual tmp<volScalarField> hf() const = 0;
};


/*---------------------------------------------------------------------------*\
           Class ubPsiMulticomponentThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class ubPsiMulticomponentThermo::implementation
:
    virtual public ubPsiMulticomponentThermo
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
        inline DimensionedFieldListAndSlicer<scalar, volMesh> Yslicer
        (
            const fvSource& model,
            const volScalarField::Internal& source
        ) const;

        //- Get the composition of a source cell
        inline scalarFieldListSlice sourceCellComposition
        (
            const DimensionedFieldListAndSlicer<scalar, volMesh>& Yslicer,
            const label i
        ) const;

        //- Get the slicer for the given source
        inline PtrList<scalarField> Yslicer
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Get the composition of a source cell
        inline scalarFieldListSlice sourceCellComposition
        (
            const PtrList<scalarField>& Yslicer,
            const label i
        ) const;
};


/*---------------------------------------------------------------------------*\
                Class ubPsiMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class ubPsiMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiThermo::implementation,
    public ubPsiMulticomponentThermo::implementation
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
            ubPsiMulticomponentThermo::implementation
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

#include "ubPsiMulticomponentThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
