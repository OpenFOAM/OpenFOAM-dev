/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::valueMulticomponentMixture

Description
    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to thermodynamic properties and mole-fraction weighted mixing to
    transport properties.

SourceFiles
    valueMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef valueMulticomponentMixture_H
#define valueMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class valueMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class valueMulticomponentMixture
:
    public multicomponentMixture<ThermoType>
{
public:

    // Public Classes

        //- Mixing type for thermodynamic properties
        class thermoMixtureType
        {
            // Private Data

                //- List of specie thermo
                const PtrList<ThermoType>& specieThermos_;

                //- List of mass fractions
                mutable List<scalar> Y_;

                //- Calculate a mass-fraction-weighted property
                template<class Method, class ... Args>
                scalar massWeighted
                (
                    Method psiMethod,
                    const Args& ... args
                ) const;

                //- Calculate a harmonic mass-fraction-weighted property
                template<class Method, class ... Args>
                scalar harmonicMassWeighted
                (
                    Method psiMethod,
                    const Args& ... args
                ) const;

                //- Limit the given temperature
                scalar limit(const scalar T) const;


        public:

            friend class valueMulticomponentMixture;


            // Constructors

                //- Construct from list of specie thermo
                thermoMixtureType(const PtrList<ThermoType>& specieThermos)
                :
                    specieThermos_(specieThermos),
                    Y_(specieThermos.size())
                {}


            // Fundamental properties

                //- Molecular weight [kg/kmol]
                scalar W() const;

                //- Return density [kg/m^3]
                scalar rho(scalar p, scalar T) const;

                //- Return compressibility [s^2/m^2]
                scalar psi(scalar p, scalar T) const;

                // Heat capacity at constant pressure [J/kg/K]
                scalar Cp(const scalar p, const scalar T) const;

                // Heat capacity at constant volume [J/kg/K]
                scalar Cv(const scalar p, const scalar T) const;

                // Sensible enthalpy [J/kg]
                scalar hs(const scalar p, const scalar T) const;

                // Absolute enthalpy [J/kg]
                scalar ha(const scalar p, const scalar T) const;

                // Enthalpy of formation [J/kg]
                scalar hf() const;


            // Mass specific derived properties

                //- Heat capacity at constant pressure/volume [J/kg/K]
                scalar Cpv(const scalar p, const scalar T) const;

                //- Gamma = Cp/Cv []
                scalar gamma(const scalar p, const scalar T) const;

                //- Enthalpy/Internal energy [J/kg]
                scalar he(const scalar p, const scalar T) const;


            // Energy->temperature  inversion functions

                //- Temperature from enthalpy or internal energy
                //  given an initial temperature T0
                scalar The
                (
                    const scalar he,
                    const scalar p,
                    const scalar T0
                ) const;
        };

        //- Mixing type for transport properties
        class transportMixtureType
        {
            // Private Data

                //- List of specie thermo
                const PtrList<ThermoType>& specieThermos_;

                //- List of mole fractions
                mutable List<scalar> X_;

                //- Calculate a mole-fraction-weighted property
                template<class Method, class ... Args>
                scalar moleWeighted
                (
                    Method psiMethod,
                    const Args& ... args
                ) const;


        public:

            friend class valueMulticomponentMixture;


            // Constructors

                //- Construct from list of specie thermo
                transportMixtureType(const PtrList<ThermoType>& specieThermos)
                :
                    specieThermos_(specieThermos),
                    X_(specieThermos.size())
                {}


            // Transport properties

                //- Dynamic viscosity [kg/m/s]
                scalar mu(const scalar p, const scalar T) const;

                //- Thermal conductivity [W/m/K]
                scalar kappa(const scalar p, const scalar T) const;
        };


private:

    // Private Data

        //- Mutable storage for the cell/face mixture thermo data
        mutable thermoMixtureType thermoMixture_;

        //- Mutable storage for the cell/face mixture transport data
        mutable transportMixtureType transportMixture_;


public:

    // Constructors

        //- Construct from a dictionary
        valueMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        valueMulticomponentMixture
        (
            const valueMulticomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~valueMulticomponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "valueMulticomponentMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valueMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
