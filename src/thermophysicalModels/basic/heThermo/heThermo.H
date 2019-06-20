/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::heThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    heThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heThermo_H
#define heThermo_H

#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class MixtureType>
class heThermo
:
    public BasicThermo,
    public MixtureType
{
protected:

    // Protected data

        //- Energy field
        volScalarField he_;


    // Protected Member Functions

        // Enthalpy/Internal energy

            //- Correct the enthalpy/internal energy field boundaries
            void heBoundaryCorrection(volScalarField& he);


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        heThermo(const heThermo<BasicThermo, MixtureType>&);


        //- Initialize heThermo
        void init();


public:

    // Constructors

        //- Construct from mesh
        heThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and dictionary
        heThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Destructor
    virtual ~heThermo();


    // Member Functions

        //- Return the compostion of the mixture
        virtual typename MixtureType::basicMixtureType&
        composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const typename MixtureType::basicMixtureType&
        composition() const
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


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }


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
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity for temperature of mixture
            //  for patch [W/m/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const;

            //- Thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "heThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
