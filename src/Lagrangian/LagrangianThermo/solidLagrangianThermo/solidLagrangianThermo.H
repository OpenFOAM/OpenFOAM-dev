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
    Foam::solidLagrangianThermo

Description
    Base-class for solid Lagrangian thermodynamic models

SourceFiles
    solidLagrangianThermo.C
    solidLagrangianThermoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidLagrangianThermo_H
#define solidLagrangianThermo_H

#include "SolidLagrangianThermo.H"
#include "pureLagrangianThermo.H"
#include "uniformDimensionedFields.H"
#include "uniformGeometricFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solidLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class solidLagrangianThermo
:
    virtual public basicLagrangianThermo
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
            SolidLagrangianThermo
            <
                BasicLagrangianThermo
                <
                    MixtureType,
                    composite
                >
            >;

        //- The derived name
        static word derivedThermoName()
        {
            return "heSolidThermo";
        }


    //- Runtime type information
    TypeName("solidLagrangianThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidLagrangianThermo,
        LagrangianMesh,
        (const LagrangianMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Select a solid thermo
        static autoPtr<solidLagrangianThermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~solidLagrangianThermo();


    // Member Functions

        //- Initialise state
        virtual void initialise();

        //- Update the pressure (does nothing)
        virtual void correctPressure(const LagrangianSubMesh&);
};


/*---------------------------------------------------------------------------*\
             Class solidLagrangianThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class solidLagrangianThermo::implementation
:
    virtual public solidLagrangianThermo
{
protected:

    // Protected Member Data

        //- Pressure [Pa]
        //  Note: This value should never be used. Solid thermo should only be
        //  instantiated with thermo models that do not depend on pressure.
        //  This uniform field takes a value of NaN, so that if any thermo
        //  models that do depend on pressure are used then the code will exit.
        uniformGeometricScalarField p_;


    // Protected Member Functions

        //- Pressure for a sub-mesh [Pa]
        tmp<uniformDimensionedScalarField> p
        (
            const LagrangianSubMesh&
        ) const;

        //- Pressure for an injection [Pa]
        tmp<uniformDimensionedScalarField> p
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const LagrangianMesh&, const word&);


    //- Destructor
    virtual ~implementation();
};


/*---------------------------------------------------------------------------*\
             Class solidLagrangianThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class solidLagrangianThermo::composite
:
    public basicLagrangianThermo::implementation,
    public pureLagrangianThermo,
    public solidLagrangianThermo::implementation
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        composite
        (
            const dictionary& dict,
            const LagrangianMesh& mesh,
            const word& phaseName
        )
        :
            basicLagrangianThermo::implementation(dict, mesh, phaseName),
            solidLagrangianThermo::implementation(dict, mesh, phaseName)
        {}

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
            composite(dict, mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
