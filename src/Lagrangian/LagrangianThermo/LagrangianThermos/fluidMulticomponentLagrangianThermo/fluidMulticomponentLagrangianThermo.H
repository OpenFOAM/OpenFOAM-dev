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
    Foam::fluidMulticomponentLagrangianThermo

Description
    Base class for multicomponent fluid Lagrangian thermodynamic models

SourceFiles
    fluidMulticomponentLagrangianThermo.C
    fluidMulticomponentLagrangianThermos.C

\*---------------------------------------------------------------------------*/

#ifndef fluidMulticomponentLagrangianThermo_H
#define fluidMulticomponentLagrangianThermo_H

#include "fluidLagrangianThermo.H"
#include "multicomponentLagrangianThermo.H"
#include "MulticomponentLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class fluidMulticomponentLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentLagrangianThermo
:
    virtual public fluidLagrangianThermo,
    virtual public multicomponentLagrangianThermo
{
public:

    // Public Classes

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            MulticomponentLagrangianThermo
            <
                FluidLagrangianThermo
                <
                    BasicLagrangianThermo
                    <
                        MixtureType,
                        composite
                    >
                >
            >;

        //- The derived name
        static word derivedThermoName()
        {
            return "heRhoThermo";
        }


    //- Runtime type information
    TypeName("fluidMulticomponentLagrangianThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidMulticomponentLagrangianThermo,
        LagrangianMesh,
        (const LagrangianMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Select a fluid multicomponent thermo
        static autoPtr<fluidMulticomponentLagrangianThermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~fluidMulticomponentLagrangianThermo();
};


/*---------------------------------------------------------------------------*\
          Class fluidMulticomponentLagrangianThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class fluidMulticomponentLagrangianThermo::composite
:
    public basicLagrangianThermo::implementation,
    public fluidLagrangianThermo::implementation,
    public multicomponentLagrangianThermo::implementation,
    public fluidMulticomponentLagrangianThermo
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const LagrangianMesh& mesh,
            const word& phaseName
        )
        :
            basicLagrangianThermo::implementation(dict, mesh, phaseName),
            fluidLagrangianThermo::implementation(dict, mesh, phaseName, T_),
            multicomponentLagrangianThermo::implementation
            (
                dict,
                mixture.species(),
                mesh,
                phaseName
            )
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
