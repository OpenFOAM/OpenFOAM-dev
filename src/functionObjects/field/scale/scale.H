/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::functionObjects::scale

Description
    Multiplies a field by a scaling factor.

    The operation can be applied to any volume or surface fields generating a
    volume or surface scalar field.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    scale.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scale_H
#define functionObjects_scale_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class scale Declaration
\*---------------------------------------------------------------------------*/

class scale
:
    public fieldExpression
{
    // Private Data

        //- Scale factor
        scalar scale_;


    // Private Member Functions

        //- Calculate the scale of the field and register the result
        template<class Type>
        bool calcScale();

        //- Calculate the scale of the field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("scale");


    // Constructors

        //- Construct from Time and dictionary
        scale
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~scale();


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "scaleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
