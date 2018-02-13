/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::chemistryReductionMethods::DAC

Description
    The Dynamic Adaptive Chemistry (DAC) method [1] simplify the chemistry
    using the matrix rAB defined by (DRGEP algorithm [2])

            |sum_i=1->Nr vAi wi dBi|
    rAB = --------------------------- ,
                max(PA, CA)

    PA = sum_i=1->Nr (max (0, vAi wi))  -> production of species A

    CA = sum_i=1->Nr (max (0, -vAi wi)) -> consumption of species A

    where i is the reaction index, Nr the number of reactions, vAi is the net
    stoechiometric coefficient of species A in the ith reaction (vAi = v''-v')
    , wi is the progress variable of reaction i and dBi equals 1 if reaction i
    involves B and O otherwise.
    rAB show the error introduced to the production rates of A when B and all
    the reactions including it are removed. It is computed as in [3] so that
    the algorithm is O(Nr).

    DAC uses a initial set of species that represents the major parts of the
    combustion mechanism, i.e. H2/O2, fuel decomposition and CO2 production.
    Usually, it includes the fuel, HO2 and CO. Then it computes the dependence
    of these set to the other species. This is done by introducing R-value
    defined by

    R_V0 (V) = max_SP(product(rij)),

    where SP is the set of all possible paths leading from V0 to V and
    product(rij) is the chain product of the weights of the edges along the
    given path. The R-value for the initial set species is 1.

    When the R-value of a species is larger than a user-defined tolerance
    then the species is included in the simplified mechanism. Otherwise,
    the species is removed along with all the reactions including it.

    During this process, instead of looking over all species like described
    in [1], the algorithm implemented here creates dynamic list to retain
    the initialized edges only (see [3]).

    References:
    \verbatim
        [1] Liang, L., Stevens, J. G., & Farrell, J. T. (2009).
        A dynamic adaptive chemistry scheme for reactive flow computations.
        Proceedings of the Combustion Institute, 32(1), 527-534.

        [2] Pepiot-Desjardins, P., & Pitsch, H. (2008).
        An efficient error-propagation-based reduction method for large
        chemical kinetic mechanisms.
        Combustion and Flame, 154(1), 67-81.

        [3] Lu, T., & Law, C. K. (2006).
        Linear time reduction of large kinetic mechanisms with directed
        relation graph: n-Heptane and iso-octane.
        Combustion and Flame, 144(1), 24-36.
    \endverbatim

SourceFiles
    DAC.C

\*---------------------------------------------------------------------------*/

#ifndef DAC_H
#define DAC_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                           Class DAC Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class DAC
:
    public chemistryReductionMethod<CompType, ThermoType>
{
    // Private data

        //- List of label for the search initiating set
        labelList searchInitSet_;

        scalar zprime_;
        label nbCLarge_;
        List<label> sC_,sH_,sO_,sN_;
        label CO2Id_,COId_,HO2Id_,H2OId_,NOId_;
        Switch automaticSIS_;
        scalar phiTol_;
        scalar NOxThreshold_;
        wordList fuelSpecies_;
        scalarField fuelSpeciesProp_;
        List<label> fuelSpeciesID_;
        word CO2Name_, COName_, HO2Name_, H2OName_, NOName_;
        Switch forceFuelInclusion_;

public:

    //- Runtime type information
    TypeName("DAC");


    // Constructors

        //- Construct from components
        DAC
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    //- Destructor
    virtual ~DAC();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DAC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
