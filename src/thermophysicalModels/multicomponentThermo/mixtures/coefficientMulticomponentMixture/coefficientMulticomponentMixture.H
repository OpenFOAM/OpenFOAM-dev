/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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

        //- Construct from dictionary, mesh and phase name
        coefficientMulticomponentMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

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

        const thermoMixtureType& cellThermoMixture(const label celli) const;

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label patchi,
            const label facei
        ) const;

        const transportMixtureType& cellTransportMixture
        (
            const label celli
        ) const
        {
            return cellThermoMixture(celli);
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return patchFaceThermoMixture(patchi, facei);
        }

        const transportMixtureType& cellTransportMixture
        (
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }
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
