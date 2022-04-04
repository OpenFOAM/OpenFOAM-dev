/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::lagrangianDistributionMap

Description
    Class containing mesh-to-mesh mapping information for particles

SourceFiles
    lagrangianDistributionMap.C

\*---------------------------------------------------------------------------*/

#ifndef lagrangianDistributionMap_H
#define lagrangianDistributionMap_H

#include "distributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class lagrangianDistributionMap Declaration
\*---------------------------------------------------------------------------*/

class lagrangianDistributionMap
{
    // Private Data

        //- Map to distribute particles
        const distributionMap particleMap_;

        //- Per element in subsetted mesh the cell label
        const labelListList constructCellLabels_;


public:

    // Constructors

        //- Construct from components
        lagrangianDistributionMap
        (
            const label nNewParticles,
            labelListList&& subParticleMap,
            labelListList&& constructParticleMap,
            labelListList&& constructCellLabels
        )
        :
            particleMap_(nNewParticles, subParticleMap, constructParticleMap),
            constructCellLabels_(constructCellLabels)
        {}


    // Member Functions

        // Access

            //- Distribution map
            const distributionMap& particleMap() const
            {
                return particleMap_;
            }

            //- Per received particle the destination cell label
            const labelListList& constructCellLabels() const
            {
                return constructCellLabels_;
            }


        // Edit

            //- Distribute list of lagrangian data
            template<class T>
            void distributeLagrangianData(List<T>& lst) const
            {
                particleMap_.distribute(lst);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
