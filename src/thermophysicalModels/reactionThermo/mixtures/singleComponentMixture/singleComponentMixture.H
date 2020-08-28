/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::singleComponentMixture

Description
    Single component mixture

SourceFiles
    singleComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef singleComponentMixture_H
#define singleComponentMixture_H

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class singleComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class singleComponentMixture
:
    public basicSpecieMixture
{

public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;

    //- Mixing type for thermodynamic properties
    typedef ThermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef ThermoType transportMixtureType;


private:

    // Private Data

        //- Thermo model
        ThermoType mixture_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        singleComponentMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~singleComponentMixture();


    // Member Functions

        //- Return the instantiated type name
        //  Note: This defines the same name as pureMixture so that
        //  non-reacting and reacting solvers are compatible
        static word typeName()
        {
            return "pureMixture<" + ThermoType::typeName() + '>';
        }

        const thermoMixtureType& cellThermoMixture(const label) const
        {
            return mixture_;
        }

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& cellTransportMixture(const label) const
        {
            return mixture_;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& cellTransportMixture
        (
            const label,
            const thermoMixtureType&
        ) const
        {
            return mixture_;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label,
            const thermoMixtureType&
        ) const
        {
            return mixture_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& specieThermo(const label speciei) const
        {
            #ifdef FULLDEBUG
            if (speciei != 0)
            {
                FatalErrorInFunction
                    << "Specie index must be zero for a single component "
                    << "mixture" << exit(FatalError);
            }
            #endif
            return mixture_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
