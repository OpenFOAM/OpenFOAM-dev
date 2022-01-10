/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::basicThermo

Description
    Base-class for fluid and solid thermodynamic properties.

    The basicThermo class is the pure virtual interface, plus static data and
    functions, including the selection table. It is default-constructable and
    holds no non-static data.

    The basicThermo::implementation class contains the class data, non-default
    construction and other implementation details.

    This structure allows for arbitrary mixing of the interfaces. Derived
    levels should implement a similar separation of interface and
    implementation. All interface classes should be default-constructable, and
    should be inherited virtually into whatever interfaces are required.
    Implementation classes should virtually inherit from their corresponding
    interface class and *not* a lower level implementation class.

    In derived levels that are complete enough to act as the base of a
    thermodynamic instantiation, an additional *::composite sub class should be
    defined, which (non-virtually) inherits from all the *::implementation
    classes that it needs. This way, the composite class composes an interface
    for which diamond patterns are resolved by virtual inheritance, whilst also
    obtaining a single copy of all the implementation classes it needs to form
    a complete implementation. The use of virtual inheritance does not result
    in additional constructor calls propagating further down the hierarchy
    (into heThermo and similar) because all virtually inherited interface
    classes are default constructable.

SourceFiles
    basicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "physicalProperties.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicThermo Declaration
\*---------------------------------------------------------------------------*/

class basicThermo
{
protected:

    // Protected Static Member Functions

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Get the constructor iterator for the given thermo dictionary and
        //  entries and constructor iterator table
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoTypeDict,
            Table* tablePtr,
            const int nCmpt,
            const char* cmptNames[],
            const word& thermoTypeName
        );

        //- Get the constructor iterator for the given thermo dictionary and
        //  constructor iterator table
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );


    // Protected Member Functions

        //- Enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryTypes();

        //- Enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary types
        wordList heBoundaryBaseTypes();


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("basicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Static Member Functions

        //- Name of a property for a given phase
        static word phasePropertyName(const word& name, const word& phaseName)
        {
            return IOobject::groupName(name, phaseName);
        }

        //- Lookup the thermo associated with the given patch field
        static const basicThermo& lookupThermo(const fvPatchScalarField& pf);

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const int nCmpt
        );

        //- Split name of thermo package into a list of named components names
        static List<Pair<word>> thermoNameComponents
        (
            const word& thermoName
        );


    // Selectors

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        //- Specialisation of the Generic New for basicThermo
        static autoPtr<basicThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~basicThermo();


    // Member Functions

        //- Properties dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties() = 0;

        //- Phase name
        virtual const word& phaseName() const = 0;

        //- Name of a property for a given phase
        word phasePropertyName(const word& name) const
        {
            return phasePropertyName(name, phaseName());
        }

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&
        ) const;

        //- Update properties
        virtual void correct() = 0;

        //- Name of the thermo physics
        virtual word thermoName() const = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const = 0;


        // Access to thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Old-time density [kg/m^3]
            virtual tmp<volScalarField> rho0() const = 0;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he() = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const = 0;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const = 0;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> ha() const = 0;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Absolute enthalpy for patch [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const = 0;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const = 0;


        // Access to transport state variables

            //- Thermal conductivity of mixture [W/m/K]
            virtual const volScalarField& kappa() const = 0;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const = 0;

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& T() = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;


        // Fields derived from transport state variables

            //- Thermal diffusivity of energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const = 0;

            //- Thermal diffusivity of energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class basicThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class basicThermo::implementation
:
    virtual public basicThermo,
    public physicalProperties
{
protected:

    // Protected data

        //- Phase-name
        const word& phaseName_;


        // Fields

            //- Temperature [K]
            volScalarField T_;

            //- Thermal conductivity [W/m/K]
            volScalarField kappa_;


        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt_;


public:

    // Static Member data

        //- Name of the thermophysical properties dictionary
        using physicalProperties::typeName;


    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;



    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Properties dictionary
        virtual const IOdictionary& properties() const
        {
            return *this;
        }

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties()
        {
            return *this;
        }

        //- Phase name
        virtual const word& phaseName() const
        {
            return phaseName_;
        }

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return dpdt_;
        }


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& T();


        // Access to transport state variables

            //- Thermal conductivity of mixture [W/m/K]
            virtual const volScalarField& kappa() const;


        //- Read thermophysical properties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
