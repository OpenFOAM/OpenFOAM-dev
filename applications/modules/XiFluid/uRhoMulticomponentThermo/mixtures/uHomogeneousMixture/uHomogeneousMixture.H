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
    Foam::uHomogeneousMixture

Description
    Unburnt gas mixture for homogeneous lean, stoichiometric or rich conditions

SourceFiles
    uHomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uHomogeneousMixture_H
#define uHomogeneousMixture_H

#include "dimensionedTypes.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class uHomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class uHomogeneousMixture
{
public:

    // Public Typedefs

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


private:

    // Private Data

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi_;

        //- Reactant thermodynamic model
        thermoType reactants_;


public:

    // Constructors

        //- Construct from a dictionary
        uHomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        uHomogeneousMixture(const uHomogeneousMixture<ThermoType>&) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "uHomogeneousMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the table of species
        const speciesTable& species() const
        {
            return speciesTable::null;
        }

        //- Access the specie active flags
        const boolList& speciesActive() const
        {
            NotImplemented;
            return boolList::null();
        }

        //- Access the specie active flags
        boolList& speciesActive()
        {
            NotImplemented;
            return const_cast<boolList&>(boolList::null());
        }

        inline const ThermoType& specieThermo(const label speciei) const
        {
            return reactants_;
        }

        //- Return the fuel-oxidant equivalence ratio
        scalar Phi(const scalarFieldListSlice& Yu) const
        {
            return Phi_;
        }

        PtrList<volScalarField::Internal> prompt
        (
            const PtrList<volScalarField>& Yu
        ) const
        {
            NotImplemented;
            return PtrList<volScalarField::Internal>();
        }

        void reset
        (
            const volScalarField& b,
            PtrList<volScalarField>& Yu,
            const volScalarField& c,
            const PtrList<volScalarField>& Yb
        ) const
        {
            FatalErrorInFunction
                << "EGR not supported by " << typeName()
                << exit(FatalError);
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const
        {
            return reactants_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const
        {
            return reactants_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const
        {
            return reactants_;
        }

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uHomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "uHomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
