/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::functionObjects::shearStress

Description
    Calculates and writes the shear-stress as the volSymmTensorField field
    'shearStress'.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    shearStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_shearStress_H
#define functionObjects_shearStress_H

#include "fieldExpression.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class shearStress Declaration
\*---------------------------------------------------------------------------*/

class shearStress
:
    public fvMeshFunctionObject,
    public writeLocalObjects
{
    // Private Data

        //- The name of the phase
        word phaseName_;


public:

    //- Runtime type information
    TypeName("shearStress");


    // Constructors

        //- Construct from Time and dictionary
        shearStress
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~shearStress();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Calculate the shearStress field
        virtual bool execute();

        //- Do nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
