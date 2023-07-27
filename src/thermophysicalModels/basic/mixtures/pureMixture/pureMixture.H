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
    Foam::pureMixture

Description
    Pure mixture model. This does no mixing, it just returns the single
    underlying thermo model.

SourceFiles
    pureMixture.C

\*---------------------------------------------------------------------------*/

#ifndef pureMixture_H
#define pureMixture_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pureMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class pureMixture
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

        //- Thermo model
        ThermoType mixture_;


public:

    // Constructors

        //- Construct from a dictionary
        pureMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        pureMixture(const pureMixture<ThermoType>&) = delete;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "pureMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture(const nil) const
        {
            return mixture_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture(const nil) const
        {
            return mixture_;
        }

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const nil,
            const thermoMixtureType&
        ) const
        {
            return mixture_;
        }

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pureMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "pureMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
