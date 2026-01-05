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
    Foam::basicLagrangianThermo

Description
    Base-class for Lagrangian fluid and solid thermodynamic models.

    The design of Lagrangian thermodynamic models is the same as that for
    finite-volume. See basicThermo for an explanation of how the interface
    class, and the implementation and composite sub-classes fit together.

See also
    Foam::basicThermo

SourceFiles
    basicLagrangianThermo.C
    basicLagrangianThermoTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef basicLagrangianThermo_H
#define basicLagrangianThermo_H

#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class basicLagrangianThermo Declaration
\*---------------------------------------------------------------------------*/

class basicLagrangianThermo
{
protected:

    // Protected Member Functions

        //- Internal energy field boundary types
        wordList eBoundaryTypes() const;

        //- Internal energy field boundary base types
        wordList eBoundaryBaseTypes() const;

        //- Field sources types
        template<class LagrangianFieldSourceType>
        static HashTable<word> sourcesTypes
        (
            const LagrangianScalarDynamicField& T
        );


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("basicLagrangianThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicLagrangianThermo,
        LagrangianMesh,
        (const LagrangianMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Select thermo of a given derived type
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );

        //- Select a basic thermo
        static autoPtr<basicLagrangianThermo> New
        (
            const LagrangianMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~basicLagrangianThermo();


    // Member Functions

        //- Access the properties dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Modify the properties dictionary
        virtual IOdictionary& properties() = 0;

        //- Access the mesh
        virtual const LagrangianMesh& mesh() const = 0;

        //- Access the phase name
        virtual const word& phaseName() const = 0;

        //- Initialise state
        virtual void initialise() = 0;

        //- Update the pressure
        virtual void correctPressure(const LagrangianSubMesh&) = 0;

        //- Update state
        virtual void correct(const LagrangianSubMesh&) = 0;

        //- Return the full name of the thermodynamic model
        virtual word thermoName() const = 0;


        // Molecular properties

            //- Molecular weight [kg/kmol]
            virtual tmp<LagrangianSubScalarField> W
            (
                const LagrangianSubMesh&
            ) const = 0;


        // Thermodynamic State

            //- Temperature [K]
            virtual const LagrangianScalarDynamicField& T() const = 0;

            //- Modify the Temperature [K]
            virtual LagrangianScalarDynamicField& T() = 0;

            //- Density [kg/m^3]
            virtual const LagrangianScalarDynamicField& rho() const = 0;

            //- Modify the Density [kg/m^3]
            virtual LagrangianScalarDynamicField& rho() = 0;

            //- Internal energy [J/kg]
            virtual const LagrangianScalarDynamicField& e() const = 0;

            //- Modify the internal energy [J/kg]
            virtual LagrangianScalarDynamicField& e() = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const LagrangianScalarDynamicField& Cv() const = 0;


        // Derived Thermodynamic Properties

            //- Density for an injection [kg/m^3]
            virtual tmp<LagrangianSubScalarField> rho
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;

            //- Internal energy for an injection [J/kg]
            virtual tmp<LagrangianSubScalarField> e
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;

            //- Heat capacity at constant volume for an injection [J/kg/K]
            virtual tmp<LagrangianSubScalarField> Cv
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;


        // Transport State

            //- Thermal conductivity [W/m/K]
            virtual const LagrangianScalarDynamicField& kappa() const = 0;


        // Derived Transport Properties

            //- Thermal conductivity for an injection [W/m/K]
            virtual tmp<LagrangianSubScalarField> kappa
            (
                const LagrangianSubScalarField& T,
                const LagrangianInjection&
            ) const = 0;
};


/*---------------------------------------------------------------------------*\
             Class basicLagrangianThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class basicLagrangianThermo::implementation
:
    virtual public basicLagrangianThermo
{
    // Private Member Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Reference to the phase name
        const word& phaseName_;


protected:

    // Protected Member Data

        //- Temperature [K]
        LagrangianScalarDynamicField T_;

        //- Density [kg/m^3]
        LagrangianScalarDynamicField rho_;

        //- Heat capacity at constant volume [J/kg/K]
        LagrangianScalarDynamicField Cv_;

        //- Thermal conductivity [W/m/K]
        LagrangianScalarDynamicField kappa_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const LagrangianMesh&, const word&);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Return const access to the mesh
        virtual const LagrangianMesh& mesh() const;

        //- Phase name
        virtual const word& phaseName() const;


        // Thermodynamic State

            //- Temperature [K]
            virtual const LagrangianScalarDynamicField& T() const;

            //- Modify the Temperature [K]
            virtual LagrangianScalarDynamicField& T();

            //- Density [kg/m^3]
            virtual const LagrangianScalarDynamicField& rho() const;

            //- Modify the Density [kg/m^3]
            virtual LagrangianScalarDynamicField& rho();

            //- Heat capacity at constant volume [J/kg/K]
            virtual const LagrangianScalarDynamicField& Cv() const;


        // Transport State

            //- Thermal conductivity [W/m/K]
            virtual const LagrangianScalarDynamicField& kappa() const;


        //- Read properties dictionary
        virtual void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicLagrangianThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
