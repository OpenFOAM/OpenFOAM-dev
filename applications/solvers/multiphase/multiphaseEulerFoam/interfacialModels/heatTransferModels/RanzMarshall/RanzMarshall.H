/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::heatTransferModels::RanzMarshall

Description

SourceFiles
    RanzMarshall.C

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshall_H
#define RanzMarshall_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                           Class RanzMarshall Declaration
\*---------------------------------------------------------------------------*/

class RanzMarshall
:
    public heatTransferModel
{

public:

    //- Runtime type information
    TypeName("RanzMarshall");


    // Constructors

        //- Construct from components
        RanzMarshall
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~RanzMarshall();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ur) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
