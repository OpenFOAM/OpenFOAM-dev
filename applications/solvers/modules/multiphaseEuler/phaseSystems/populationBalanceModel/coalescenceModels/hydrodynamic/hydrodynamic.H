/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::hydrodynamic

Description
    Hydrodynamic kernel.

SourceFiles
    hydrodynamic.C

\*---------------------------------------------------------------------------*/

#ifndef hydrodynamic_H
#define hydrodynamic_H

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                        Class hydrodynamic Declaration
\*---------------------------------------------------------------------------*/

class hydrodynamic
:
    public coalescenceModel
{
public:

    //- Runtime type information
    TypeName("hydrodynamic");

    // Constructor

        hydrodynamic
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~hydrodynamic()
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
