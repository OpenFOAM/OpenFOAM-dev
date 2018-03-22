/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::ReactingPhaseModel

Description
    Class which represents phases with volumetric reactions. Returns the
    reaction rate and heat.

SourceFiles
    ReactingPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingPhaseModel_H
#define ReactingPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ReactingPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class ReactionType>
class ReactingPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Reaction model
        autoPtr<ReactionType> reaction_;


public:

    // Constructors

        ReactingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~ReactingPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Return the species fraction equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const;

        //- Return heat release rate
        virtual tmp<volScalarField> Qdot() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
