/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::valueMultiComponentMixture

Description
    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to thermodynamic properties and mole-fraction weighted mixing to
    transport properties.

SourceFiles
    valueMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef valueMultiComponentMixture_H
#define valueMultiComponentMixture_H

#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class valueMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class valueMultiComponentMixture
:
    public multiComponentMixture<ThermoType>
{

public:

    class thermoMixture
    {
        //- List of specie thermo
        const PtrList<ThermoType>& specieThermos_;

        //- List of mass fractions
        mutable List<scalar> Y_;

        template<class Method, class ... Args>
        scalar massWeighted(Method psiMethod, const Args& ... args) const;

        template<class Method, class ... Args>
        scalar harmonicMassWeighted
        (
            Method psiMethod,
            const Args& ... args
        ) const;

        scalar limit(const scalar T) const;


    public:

        friend class valueMultiComponentMixture;

        // Constructors

            thermoMixture
            (
                const PtrList<ThermoType>& specieThermos
            )
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
            scalar Hs(const scalar p, const scalar T) const;

            // Absolute enthalpy [J/kg]
            scalar Ha(const scalar p, const scalar T) const;

            // Enthalpy of formation [J/kg]
            scalar Hf() const;


        // Mass specific derived properties

            //- Heat capacity at constant pressure/volume [J/kg/K]
            scalar Cpv(const scalar p, const scalar T) const;

            //- Gamma = Cp/Cv []
            scalar gamma(const scalar p, const scalar T) const;

            //- Ratio of heat capacity at constant pressure to that at
            //  constant pressure/volume []
            scalar CpByCpv(const scalar p, const scalar T) const;

            //- Enthalpy/Internal energy [J/kg]
            scalar HE(const scalar p, const scalar T) const;


        // Energy->temperature  inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            scalar THE
            (
                const scalar he,
                const scalar p,
                const scalar T0
            ) const;
    };


    class transportMixture
    {
        //- List of specie thermo
        const PtrList<ThermoType>& specieThermos_;

        //- List of mole fractions
        mutable List<scalar> X_;

        template<class Method, class ... Args>
        scalar moleWeighted(Method psiMethod, const Args& ... args) const;


    public:

        friend class valueMultiComponentMixture;

        transportMixture
        (
            const PtrList<ThermoType>& specieThermos
        )
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


    //- Mixing type for thermodynamic properties
    typedef thermoMixture thermoMixtureType;

    //- Mixing type for transport properties
    typedef transportMixture transportMixtureType;


private:

    // Private Data

        //- Mutable storage for the cell/face mixture thermo data
        mutable thermoMixtureType thermoMixture_;

        //- Mutable storage for the cell/face mixture transport data
        mutable transportMixtureType transportMixture_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        valueMultiComponentMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        valueMultiComponentMixture
        (
            const valueMultiComponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~valueMultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "valueMultiComponentMixture<" + ThermoType::typeName() + '>';
        }

        const thermoMixtureType& cellThermoMixture(const label celli) const;

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label patchi,
            const label facei
        ) const;

        const transportMixtureType& cellTransportMixture
        (
            const label celli
        ) const;

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei
        ) const;

        const transportMixtureType& cellTransportMixture
        (
            const label celli,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return cellTransportMixture(celli);
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return patchFaceTransportMixture(patchi, facei);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "valueMultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
