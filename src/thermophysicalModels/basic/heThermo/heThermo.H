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
    Foam::heThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    heThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heThermo_H
#define heThermo_H

#include "basicMixture.H"
#include "volFields.H"
#include "uniformGeometricFields.H"

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

        //- Heat capacity at constant pressure field [J/kg/K]
        volScalarField Cp_;

        // Heat capacity at constant volume field [J/kg/K]
        volScalarField Cv_;


    // Protected Member Functions

        //- Return a volScalarField of the given property
        template
        <
            class CellMixture,
            class PatchFaceMixture,
            class Method,
            class ... Args
        >
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            CellMixture cellMixture,
            PatchFaceMixture patchFaceMixture,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class CellMixture, class Method, class ... Args>
        tmp<scalarField> cellSetProperty
        (
            CellMixture cellMixture,
            Method psiMethod,
            const labelList& cells,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class PatchFaceMixture, class Method, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            PatchFaceMixture patchFaceMixture,
            Method psiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return an indirect list of a field for the given set of cells
        static UIndirectList<scalar> cellSetScalarList
        (
            const volScalarField& psi,
            const labelList& cells
        );

        //- Return an indirect list of a field for the given set of cells
        static UniformField<scalar> cellSetScalarList
        (
            const uniformGeometricScalarField& psi,
            const labelList&
        );

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);


public:

    // Constructors

        //- Construct from mesh
        heThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        heThermo(const heThermo<BasicThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~heThermo();


    // Member Functions

        //- Return the composition of the mixture
        virtual typename MixtureType::basicMixtureType&
        composition()
        {
            return *this;
        }

        //- Return the composition of the mixture
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

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
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

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
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

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;


        // Fields derived from transport state variables

            //- Thermal diffusivity of energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity of energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;


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
