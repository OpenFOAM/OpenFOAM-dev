/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Aggregation model of Dahneke (1983). Utilizes physical, i.e. collisional
    diameters.

    References:
    \verbatim
        Dahneke, B. (1983).
        Simple kinetic theory of Brownian diffusion in vapors and aerosols.
        In Theory of dispersed multiphase flow (pp. 97-133). Academic Press.
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

class brownianCollisions;
class ballisticCollisions;

/*---------------------------------------------------------------------------*\
                    Class DahnekeInterpolation Declaration
\*---------------------------------------------------------------------------*/

class DahnekeInterpolation
:
    public coalescenceModel
{
    // Private Data

        //- Model for coagulation due to brownian collisions
        autoPtr<brownianCollisions> brownian_;

        //- Rate for coagulation due to brownian collisions
        volScalarField brownianCollisionRate_;

        //- Model for coagulation due to ballistic collisions
        autoPtr<ballisticCollisions> ballistic_;

        //- Rate for coagulation due to ballistic collisions
        volScalarField ballisticCollisionRate_;


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
