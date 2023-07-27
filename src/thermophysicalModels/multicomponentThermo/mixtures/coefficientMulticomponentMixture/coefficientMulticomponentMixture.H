/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::coefficientMulticomponentMixture

Description
    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to the thermodynamic and transport coefficients.

SourceFiles
    coefficientMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef coefficientMulticomponentMixture_H
#define coefficientMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class coefficientMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class coefficientMulticomponentMixture
:
    public multicomponentMixture<ThermoType>
{
public:

    // Public Typedefs

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoMixtureType mixture_;


public:

    // Constructors

        //- Construct from a dictionary
        coefficientMulticomponentMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        coefficientMulticomponentMixture
        (
            const coefficientMulticomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~coefficientMulticomponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "multicomponentMixture<" + ThermoType::typeName() + '>';
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
    #include "coefficientMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
