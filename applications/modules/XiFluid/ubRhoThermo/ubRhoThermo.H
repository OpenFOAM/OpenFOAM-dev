/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::ubRhoThermo

Description

SourceFiles
    ubRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ubRhoThermo_H
#define ubRhoThermo_H

#include "uRhoMulticomponentThermo.H"
#include "bRhoMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ubRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class ubRhoThermo
:
    virtual public rhoFluidThermo
{
    // Private Data

        const word unburntPhaseName_;
        const word burntPhaseName_;

        volScalarField b_;
        volScalarField c_;

        autoPtr<uRhoMulticomponentThermo> uThermo_;
        autoPtr<bRhoMulticomponentThermo> bThermo_;

        volScalarField rho_;
        volScalarField psi_;
        volScalarField mu_;
        volScalarField kappa_;

        volScalarField alphau_;
        volScalarField alphab_;


public:

    //- Runtime type information
    TypeName("ubRhoThermo");


    // Constructors

        //- Construct from mesh
        ubRhoThermo(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        ubRhoThermo(const ubRhoThermo&) = delete;


    //- Destructor
    virtual ~ubRhoThermo();


    // Member Functions

        //- Regress variable
        const volScalarField& b() const
        {
            return b_;
        }

        //- Regress variable
        //  Non-const access allowed for transport equations
        volScalarField& b()
        {
            return b_;
        }

        //- Progress variable
        const volScalarField& c() const
        {
            return c_;
        }

        //- Progress variable
        //  Non-const access allowed for transport equations
        volScalarField& c()
        {
            return c_;
        }

        //- Unburnt gas volume fraction
        const volScalarField& alphau() const
        {
            return alphau_;
        }

        //- Burnt gas volume fraction
        const volScalarField& alphab() const
        {
            return alphab_;
        }

        const uRhoMulticomponentThermo& uThermo() const
        {
            return uThermo_();
        }

        uRhoMulticomponentThermo& uThermo()
        {
            return uThermo_();
        }

        const bRhoMulticomponentThermo& bThermo() const
        {
            return bThermo_();
        }

        bRhoMulticomponentThermo& bThermo()
        {
            return bThermo_();
        }

        //- Properties dictionary
        virtual const IOdictionary& properties() const;

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties();

        //- Return const access to the mesh
        virtual const fvMesh& mesh() const;

        //- Phase name
        virtual const word& phaseName() const;

        //- Update properties
        virtual void correct();

        //- Return the name of the mixture (not implemented)
        virtual word mixtureName() const;

        //- Name of the thermo physics (not implemented)
        virtual word thermoName() const;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return false;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return false;
        }

        //- Should the dpdt term be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return true;
        }

        //- Reset the mixture to an unburnt state and update EGR
        virtual void reset();


        // Molecular properties

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;


        // Thermodynamic state

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& T();

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he();

            //- Heat capacity at constant pressure [J/kg/K]
            virtual const volScalarField& Cp() const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual const volScalarField& Cv() const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual const volScalarField& Cpv() const;


        // Derived Thermodynamic Properties

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

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
            virtual tmp<volScalarField::Internal> he
            (
                const volScalarField::Internal& T,
                const fvSource& model,
                const volScalarField::Internal& source
            ) const;

            //- Enthalpy/Internal energy for source [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const fvSource& model,
                const scalarField& source,
                const labelUList& cells
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

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField::Internal> hs
            (
                const volScalarField::Internal& p,
                const volScalarField::Internal& T
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

            //- Absolute enthalpy [J/kg]
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

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Absolute enthalpy for patch [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
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
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> The
            (
                const scalarField& h,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        // Transport state

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const;

            //- Thermal conductivity of mixture [W/m/K]
            virtual const volScalarField& kappa() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
