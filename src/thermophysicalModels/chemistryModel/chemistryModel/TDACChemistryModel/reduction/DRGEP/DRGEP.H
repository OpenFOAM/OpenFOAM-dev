/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::chemistryReductionMethods::DRGEP

Description
    The DRGEP algorithm [1] is based on

            |sum_i=1->Nr vAi wi dBi|
    rAB = --------------------------- ,
                max(PA, CA)

    PA = sum_i=1->Nr (max (0, vAi wi))  -> production of species A

    CA = sum_i=1->Nr (max (0, -vAi wi)) -> consumption of species A

    where i is the reaction index, Nr the number of reactions, vAi is the net
    stoichiometric coefficient of species A in the ith reaction (vAi = v''-v')
    , wi is the progress variable of reaction i and dBi equals 1 if reaction i
    involves B and O otherwise.
    rAB show the error introduced to the production rates of A when B and all
    the reactions including it are removed. It is computed as in [2] so that
    the algorithm is O(Nr).

    DAC uses a initial set of species that represents the major parts of the
    combustion mechanism, i.e. H2/O2, fuel decomposition and CO2 production.
    Usually, it includes the fuel, HO2 and CO. Then it computes the dependence
    of these set to the other species. This is done by introducing R-value
    defined by

    R_V0 (V) = max_SP(product(rij)) ,

    where SP is the set of all possible paths leading from V0 to V and
    product(rij) is the chain product of the weights of the edges along the
    given path. The R-value for the initial set species is 1.

    When the R-value of a species is larger than a user-defined tolerance
    then the species is included in the simplified mechanism. Otherwise,
    the species is removed along with all the reactions including it.

    During this process, instead of looking over all species like described
    in [1], the algorithm implemented here creates dynamic list to retain
    the initialised edges only (see [2]).

    To avoid using the target species when they are not contributing yet or
    anymore to the system, a coefficient based on the exchange of element is
    introduced:

                NTa |PT - CT|
    alphaTa = ----------------
                    Pa

    Pa = sum_speciesS NSa max(0, PS-CS)

    where 'a' refers to different elements present in the system
    (namely C, H, O and N for conventional hydrocarbon combustion),
    NTa is the number of element a in species T.
    When this coefficient alpha is below the specified threshold, the species is
    removed from the search initiating set. In the original paper from
    Pepiot-Desjardins et al.[2], this coefficient is further transformed to
    compute a global normalised scaling coefficient but here as it is
    dynamically computed, alpha is not introduced in the calculation of R.

    References:
    \verbatim
        [1] Pepiot-Desjardins, P., & Pitsch, H. (2008).
        An efficient error-propagation-based reduction method for large
        chemical kinetic mechanisms.
        Combustion and Flame, 154(1), 67-81.

        [2] Lu, T., & Law, C. K. (2006).
        Linear time reduction of large kinetic mechanisms with directed
        relation graph: n-Heptane and iso-octane.
        Combustion and Flame, 144(1), 24-36.
    \endverbatim

SourceFiles
    DRGEP.C

\*---------------------------------------------------------------------------*/

#ifndef DRGEP_H
#define DRGEP_H

#include "chemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryReductionMethods
{

/*---------------------------------------------------------------------------*\
                             Class ode Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class DRGEP
:
    public chemistryReductionMethod<ThermoType>
{
    // Private Data

        //- List of label for the search initiating set
        labelList searchInitSet_;

        List<label> sC_,sH_,sO_,sN_;
        label NGroupBased_;

public:

    //- Runtime type information
    TypeName("DRGEP");


    // Constructors

        //- Construct from components
        DRGEP
        (
            const IOdictionary& dict,
            TDACChemistryModel<ThermoType>& chemistry
        );


    //- Destructor
    virtual ~DRGEP();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DRGEP.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
