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
    Foam::bInhomogeneousMixture

Description
    Burnt gas combustion mixture for inhomogeneous lean, stoichiometric
    or rich conditions

    in which the unburnt fuel mass-fraction \c fu (mixture fraction defines
    the composition.

SourceFiles
    bInhomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef bInhomogeneousMixture_H
#define bInhomogeneousMixture_H

#include "dimensionedTypes.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bInhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class bInhomogeneousMixture
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

        //- Stoichiometric air/fuel ratio
        scalar stoicRatio_;

        //- Fuel thermodynamic model
        thermoType fuel_;

        //- Oxidant thermodynamic model
        thermoType oxidant_;

        //- Product thermodynamic model
        thermoType products_;

        //- List of specie active flags
        boolList active_;

        //- Mutable storage for the mixed thermodynamic model
        mutable thermoType mixture_;


    // Private Member Functions

        //- Return the mixture for the given composition
        const thermoType& mixture(const scalar ft) const;


public:

        //- Specie indices
        enum species { FT };


    // Constructors

        //- Construct from a dictionary
        bInhomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        bInhomogeneousMixture
        (
            const bInhomogeneousMixture<ThermoType>&
        ) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "bInhomogeneousMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the specie names
        static wordList specieNames()
        {
            return {"ft"};
        }

        //- Access the specie active flags
        const boolList& speciesActive() const
        {
            return active_;
        }

        //- Access the specie active flags
        boolList& speciesActive()
        {
            return active_;
        }

        inline const ThermoType& specieThermo(const label speciei) const
        {
            FatalErrorInFunction
                << "Cannot return specieThermo for specie " << speciei << " ft "
                << exit(FatalError);
            return products_;
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

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bInhomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "bInhomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
