/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::liftModels::noLift

Description

SourceFiles
    noLift.C

\*---------------------------------------------------------------------------*/

#ifndef noLift_H
#define noLift_H

#include "liftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                           Class noLift Declaration
\*---------------------------------------------------------------------------*/

class noLift
:
    public liftModel
{
protected:

    // Protected data

        //- Interface
        const phaseInterface interface_;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and an interface
        noLift
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~noLift();


    // Member Functions

        //- Return lift force
        virtual tmp<volVectorField> F() const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
