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
#include "multicomponentThermo.H"
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
    virtual public psiThermo,
    virtual public multicomponentThermo
{

public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType = UBPsiMulticomponentThermo
        <
            MulticomponentThermo
            <
                BasicThermo<MixtureType, composite>
            >
        >;

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

        // Mass fractions

            //- Return the fuel-oxidant equivalence ratio
            virtual tmp<volScalarField> Phi() const = 0;

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
              Class ubPsiMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class ubPsiMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public psiThermo::implementation,
    public multicomponentThermo::implementation,
    public ubPsiMulticomponentThermo
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
            multicomponentThermo::implementation
            (
                dict,
                mixture.specieNames(),
                mesh,
                phaseName,
                false
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
