/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::functionObjects::blendingFactor

Description
    Calculates and outputs the blendingFactor as used by the bended convection
    schemes.  The output is a volume field (cells) whose value is calculated via
    the maximum blending factor for any cell face.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    blendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blendingFactor_H
#define functionObjects_blendingFactor_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
:
    public fieldExpression
{
    // Private member data

        //- Name of flux field, default is "phi"
        word phiName_;


    // Private Member Functions

        //- Calculate the blending factor field
        template<class Type>
        bool calcBF();

        //- Calculate the blending factor field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        blendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blendingFactor();


    // Member Functions

        //- Read the blendingFactor data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
