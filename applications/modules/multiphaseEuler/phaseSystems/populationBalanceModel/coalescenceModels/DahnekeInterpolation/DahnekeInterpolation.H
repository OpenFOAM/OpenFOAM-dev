/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::DahnekeInterpolation

Description
    Interpolation formula of Dahneke (1983) as presented by Otto et al. (1999).
    Utilises collisional diameters.

    References:
    \verbatim
        Dahneke, B. (1983).
        Simple kinetic theory of Brownian diffusion in vapors and aerosols.
        In Theory of dispersed multiphase flow (pp. 97-133). Academic Press.
    \endverbatim

    \verbatim
        Otto, E., Fissan, H., Park, S. H., & Lee, K. W. (1999).
        The log-normal size distribution theory of Brownian aerosol coagulation
        for the entire particle size range: part II—analytical solution using
        Dahneke’s coagulation kernel.
        Journal of aerosol science, 30(1), 17-34.
    \endverbatim

SourceFiles
    DahnekeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef DahnekeInterpolation_H
#define DahnekeInterpolation_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

class BrownianCollisions;
class ballisticCollisions;

/*---------------------------------------------------------------------------*\
                    Class DahnekeInterpolation Declaration
\*---------------------------------------------------------------------------*/

class DahnekeInterpolation
:
    public coalescenceModel
{
    // Private Data

        //- Model for coagulation due to Brownian collisions
        autoPtr<BrownianCollisions> Brownian_;

        //- Rate for coagulation due to Brownian collisions
        volScalarField BrownianRate_;

        //- Model for coagulation due to ballistic collisions
        autoPtr<ballisticCollisions> ballistic_;

        //- Rate for coagulation due to ballistic collisions
        volScalarField ballisticRate_;


public:

    //- Runtime type information
    TypeName("DahnekeInterpolation");

    // Constructor

        DahnekeInterpolation
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~DahnekeInterpolation()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
