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
    Foam::diameterModels::driftModels::constantDrift

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constantDrift_H
#define constantDrift_H

#include "driftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace driftModels
{

/*---------------------------------------------------------------------------*\
                        Class constantDrift Declaration
\*---------------------------------------------------------------------------*/

class constantDrift
:
    public driftModel
{
    // Private Data

        //- Drift rate
        dimensionedScalar rate_;


public:

    //- Runtime type information
    TypeName("constant");

    // Constructor

        constantDrift
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantDrift()
    {}


    // Member Functions

        //- Add to driftRate
        virtual void addToDriftRate
        (
            volScalarField& driftRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace driftModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
