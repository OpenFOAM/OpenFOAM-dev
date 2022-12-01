/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::constSolidThermo

Description
    Uniform or non-uniform constant solid thermodynamic properties

    Each physical property can specified as either \c uniform in which case the
    value entry is read, \c zonal in which case the value entry and zone list
    are read or \c file in which case the field file in read from the constant
    directory.

Usage
    Example of uniform constant solid properties specification:
    \verbatim
        thermoType          constSolidThermo;

        rho
        {
            type        uniform;
            value       8940;
        }

        Cv
        {
            type        uniform;
            value       385;
        }

        kappa
        {
            type        uniform;
            value       380;
        }
    \endverbatim

    Example of zonal constant solid properties specification where kappa is
    different in different zones:
    \verbatim
        thermoType          constSolidThermo;

        rho
        {
            type        uniform;
            value       8940;
        }

        Cv
        {
            type        uniform;
            value       385;
        }

        kappa
        {
            type        zonal;
            value       380;

            zones
            {
                heater      560;
                insulation  100;
            }
        }
    \endverbatim

    Example of non-uniform constant solid properties specification:
    \verbatim
        thermoType          constSolidThermo;

        rho
        {
            type        file;
        }

        Cv
        {
            type        file;
        }

        kappa
        {
            type        file;
        }
    \endverbatim
    where each of the field files are read from the constant directory.

SourceFiles
    constSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class constSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class constSolidThermo
:
    public solidThermo::composite
{
    // Private data

        //- Heat capacity at constant volume [J/kg/K]
        volScalarField Cv_;

        //- Internal energy [J/kg]
        volScalarField e_;


protected:

    // Protected constructors

        //- Construct from mesh and phase name
        constSolidThermo
        (
            const fvMesh&,
            const bool readKappa,
            const word& phaseName = word::null
        );


    // Protected Member Functions

        template<class Type>
        VolField<Type> readProperty
        (
            const word& name,
            const dimensionSet& dimensions
        ) const;


public:

    //- Runtime type information
    TypeName("constSolidThermo");


    // Constructors

        //- Construct from mesh and phase name
        constSolidThermo
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~constSolidThermo();


    // Member Functions

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return type();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return true;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return true;
        }


        // Access to thermophysical state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he();

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const;


        // Access to transport state variables

            //- Return true as the thermal conductivity is isotropic
            virtual bool isotropic() const
            {
                return true;
            }

            //- Anisotropic thermal conductivity [W/m/K]
            //  Not implemented
            virtual const volVectorField& Kappa() const;


        // Fields derived from thermodynamic state variables

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
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

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Sensible enthalpy for patch [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
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

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
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


        //- Update properties
        virtual void correct();


    // I-O

        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "constSolidThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
