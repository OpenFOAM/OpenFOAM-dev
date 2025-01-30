/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::BasicThermo

Description
    Thermo implementation and storage of energy and heat capacities. Provides
    overloads of the functions defined in the basic thermo type that depend on
    the primitive thermo model.

SourceFiles
    BasicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef BasicThermo_H
#define BasicThermo_H

#include "volFields.H"
#include "physicalProperties.H"
#include "uniformGeometricFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BasicThermoName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(BasicThermo);


/*---------------------------------------------------------------------------*\
                          Class BasicThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType, class BasicThermoType>
class BasicThermo
:
    public BasicThermoName,
    public physicalProperties,
    public MixtureType,
    public BasicThermoType
{
protected:

    // Protected data

        //- Energy field
        volScalarField he_;

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        // Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;


    // Protected Member Functions

        //- Return a volScalarField of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a volScalarField::Internal of the given property
        template<class Mixture, class Method, class ... Args>
        tmp<volScalarField::Internal> volInternalScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Mixture mixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class Mixture, class Method, class ... Args>
        tmp<scalarField> cellSetProperty
        (
            Mixture mixture,
            Method psiMethod,
            const labelList& cells,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class Mixture, class Method, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            Mixture mixture,
            Method psiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property for a source
        template<class Mixture, class Method, class ... Args>
        tmp<scalarField> fieldSourceProperty
        (
            Mixture mixture,
            Method psiMethod,
            const fvSource& source,
            const Args& ... args
        ) const;

        //- Return an indirect list of a field for the given set of cells
        static UIndirectList<scalar> cellSetScalarList
        (
            const volScalarField& psi,
            const labelUList& cells
        );

        //- Return an indirect list of a field for the given set of cells
        static UniformField<scalar> cellSetScalarList
        (
            const uniformGeometricScalarField& psi,
            const labelUList&
        );

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);


public:

    // Typedefs

        //- Mixture type
        typedef MixtureType mixtureType;

        //- Basic thermo
        typedef BasicThermoType basicThermoType;


    //- Disambiguate debug switch used by derivations
    using BasicThermoName::debug;


    // Constructors

        //- Construct from mesh and phase name
        BasicThermo(const fvMesh&, const word& phaseName);

        //- Disallow default bitwise copy construction
        BasicThermo(const BasicThermo<MixtureType, BasicThermoType>&) = delete;


    //- Destructor
    virtual ~BasicThermo();


    // Member Functions

        //- Return the properties dictionary
        virtual IOdictionary& properties()
        {
            return *this;
        }

        //- Return the properties dictionary
        virtual const IOdictionary& properties() const
        {
            return *this;
        }

        //- Return the mixture
        const MixtureType& mixture() const
        {
            return *this;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return MixtureType::thermoType::typeName();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return MixtureType::thermoType::incompressible;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return MixtureType::thermoType::isochoric;
        }


        // Molecular properties

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;


        // Thermodynamic state

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const
            {
                return Cp_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const
            {
                return Cv_;
            }

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const;


        // Derived Thermodynamic Properties

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField::Internal> he
            (
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Enthalpy/Internal energy for source [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const fvSource& source
            ) const;

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField::Internal> hs
            (
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Sensible enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField::Internal> ha
            (
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
            ) const;

            //- Absolute enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const;


        // Temperature-energy inversion functions

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> The
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> The
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        //- Read thermophysical properties dictionary
        using BasicThermoType::read;

        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
