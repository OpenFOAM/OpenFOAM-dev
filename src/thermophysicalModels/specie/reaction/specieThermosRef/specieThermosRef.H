/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::specieThermosRef

Description
    Minimal wrapper around a reference to a list of specie thermos that allows
    it to passed through layers of code which do not know about (i.e., aren't
    templated on) the form of the thermo model

\*---------------------------------------------------------------------------*/

#ifndef specieThermosRef_H
#define specieThermosRef_H

#include "typeInfo.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ThermoType>
class SpecieThermosRef;

/*---------------------------------------------------------------------------*\
                          Class specieThermosRef Declaration
\*---------------------------------------------------------------------------*/

class specieThermosRef
{
public:

    //- Runtime type information
    TypeName("specieThermosRef");


    //- Destructor
    virtual ~specieThermosRef()
    {}


    // Member Functions

        //- Get the type of the species' thermos
        virtual word thermoName() const = 0;

        //- Cast to derived, given the type of the species' thermos
        template<class ThermoType>
        const SpecieThermosRef<ThermoType>& cast() const
        {
            return refCast<const SpecieThermosRef<ThermoType>>(*this);
        }
};


/*---------------------------------------------------------------------------*\
                          Class SpecieThermosRef Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class SpecieThermosRef
:
    public specieThermosRef
{
    // Private Member Data

        //- The specie thermos
        const PtrList<ThermoType>& specieThermos_;


public:

    // Constructors

        //- Construct from a reference to the specie thermos
        SpecieThermosRef(const PtrList<ThermoType>& specieThermos)
        :
            specieThermos_(specieThermos)
        {}


    // Member Functions

        //- Get the type of the species' thermos
        virtual word thermoName() const
        {
            return ThermoType::typeName();
        }

        //- Get the number of specie thermos
        label size() const
        {
            return specieThermos_.size();
        }

        //- Access a specie thermo
        const ThermoType& operator[](const label speciei) const
        {
            return specieThermos_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
