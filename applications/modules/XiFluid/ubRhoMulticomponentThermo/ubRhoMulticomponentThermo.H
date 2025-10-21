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
    Foam::ubRhoMulticomponentThermo

Description
    Base-class for combustion fluid thermodynamic properties based on
    compressibility.

See also
    Foam::basicThermo

SourceFiles
    ubRhoMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ubRhoMulticomponentThermo_H
#define ubRhoMulticomponentThermo_H

#include "rhoFluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "UBRhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ubRhoMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

class ubRhoMulticomponentThermo
:
    virtual public rhoFluidThermo,
    virtual public fluidMulticomponentThermo
{

public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType = UBRhoFluidMulticomponentThermo
        <
            FluidMulticomponentThermo
            <
                MulticomponentThermo
                <
                    BasicThermo<MixtureType, composite>
                >
            >
        >;

        //- The derived name
        static const word derivedThermoName()
        {
            return "ubRhoThermo";
        }


    //- Runtime type information
    TypeName("ubRhoMulticomponentThermo");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            ubRhoMulticomponentThermo,
            fvMesh,
            (const fvMesh& mesh, const word& phaseName),
            (mesh, phaseName)
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<ubRhoMulticomponentThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~ubRhoMulticomponentThermo();


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
              Class ubRhoMulticomponentThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class ubRhoMulticomponentThermo::composite
:
    public basicThermo::implementation,
    public rhoThermo::implementation,
    public fluidThermo::implementation,
    public multicomponentThermo::implementation,
    public ubRhoMulticomponentThermo
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
            rhoThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName),
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

#include "ubRhoMulticomponentThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
