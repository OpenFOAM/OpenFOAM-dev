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
    Foam::BHomogeneousMixture

Description
    Templated burnt gas mixture for homogeneous lean, stoichiometric or rich
    conditions

See also
    Foam::bHomogeneousMixture

SourceFiles
    BHomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef BHomogeneousMixture_H
#define BHomogeneousMixture_H

#include "bHomogeneousMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BHomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BHomogeneousMixture
:
    public bHomogeneousMixture
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

        //- Products thermodynamic model
        thermoType products_;


public:

    // Constructors

        //- Construct from a dictionary
        BHomogeneousMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        BHomogeneousMixture(const BHomogeneousMixture<ThermoType>&) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                bHomogeneousMixture::typeName
              + '<' + ThermoType::typeName() + '>';
        }

        inline const ThermoType& specieThermo(const label speciei) const
        {
            return products_;
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const
        {
            return products_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const
        {
            return products_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const
        {
            return products_;
        }

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BHomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BHomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
