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
    Foam::coefficientWilkeMulticomponentMixture

Description
    Thermophysical properties mixing class which applies mass-fraction weighted
    mixing to the thermodynamic coefficients and Wilke's equation to
    transport properties.

    Reference:
    \verbatim
        Wilke, C. R. (1950).
        A viscosity equation for gas mixtures.
        The journal of chemical physics, 18(4), 517-519.
    \endverbatim

SourceFiles
    coefficientWilkeMulticomponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef coefficientWilkeMulticomponentMixture_H
#define coefficientWilkeMulticomponentMixture_H

#include "multicomponentMixture.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coefficientWilkeMulticomponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class coefficientWilkeMulticomponentMixture
:
    public multicomponentMixture<ThermoType>
{

public:

    class transportMixture
    {
        //- List of specie thermo
        const PtrList<ThermoType>& specieThermos_;

        //- List of molecular weights
        scalarList M_;

        //- Matrix of (4/sqrt(2.0))*sqrt(1 + M_[i]/M_[j])
        scalarSquareMatrix A_;

        //- Matrix of sqrt(M_[j]/M_[i])
        scalarSquareMatrix B_;

        //- List of mole fractions
        mutable scalarList X_;

        //- List of specie viscosities
        mutable scalarList mu_;

        //- List of Wilke weights
        mutable scalarList w_;

        //- mu cache state to avoid recalculation of the Wilke weight for kappa
        mutable bool muCached_;

        //- Calculate the Wilke weights and store in w_
        void WilkeWeights(const scalar p, const scalar T) const;


    public:

        friend class coefficientWilkeMulticomponentMixture;

        // Constructors

            transportMixture
            (
                const PtrList<ThermoType>& specieThermos
            );


        // Transport properties

            //- Dynamic viscosity [kg/m/s]
            scalar mu(const scalar p, const scalar T) const;

            //- Thermal conductivity [W/m/K]
            scalar kappa(const scalar p, const scalar T) const;
    };


public:

    //- Mixing type for thermodynamic properties
    typedef typename ThermoType::thermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef transportMixture transportMixtureType;


private:

    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable thermoMixtureType mixture_;

        //- Mutable storage for the cell/face mixture transport data
        mutable transportMixtureType transportMixture_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        coefficientWilkeMulticomponentMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        coefficientWilkeMulticomponentMixture
        (
            const coefficientWilkeMulticomponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~coefficientWilkeMulticomponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "coefficientWilkeMulticomponentMixture<"
               + ThermoType::typeName() + '>';
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
        ) const;

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei,
            const thermoMixtureType& thermoMixture
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coefficientWilkeMulticomponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
