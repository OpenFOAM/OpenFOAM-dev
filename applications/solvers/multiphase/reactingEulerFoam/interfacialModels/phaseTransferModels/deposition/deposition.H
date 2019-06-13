/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::phaseTransferModels::deposition

Description
    Phase transfer model representing change from a dispersed phase to a film as
    a result of deposition onto a third phase

SourceFiles
    deposition.C

\*---------------------------------------------------------------------------*/

#ifndef deposition_H
#define deposition_H

#include "phaseTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace phaseTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class deposition Declaration
\*---------------------------------------------------------------------------*/

class deposition
:
    public phaseTransferModel
{
private:

    // Private Data

        //- The name of the phase which deposits
        const word dropletName_;

        //- The name of the phase onto which deposition occurs
        const word surfaceName_;

        //- The deposition efficiency
        const scalar efficiency_;


public:

    //- Runtime type information
    TypeName("deposition");


    // Constructors

        //- Construct from components
        deposition
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~deposition();


    // Member Functions

        //- The mass transfer rate
        virtual tmp<volScalarField> dmdt() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
