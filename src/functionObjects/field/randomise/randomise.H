/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::functionObjects::randomise

Description
    Adds a random component to a field, with a specified perturbation magnitude.

    The operation can be applied to any volume field.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    randomise.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_randomise_H
#define functionObjects_randomise_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class randomise Declaration
\*---------------------------------------------------------------------------*/

class randomise
:
    public fieldExpression
{
    // Private member data

        //- The magnitude of the purturbation
        scalar magPerturbation_;


    // Private Member Functions

        //- Calculate the randomisenitude of the field and register the result
        template<class Type>
        bool calcRandomised();

        //- Calculate the randomised field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("randomise");


    // Constructors

        //- Construct from Time and dictionary
        randomise
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~randomise();


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "randomiseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
