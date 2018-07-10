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
    Foam::functionObjects::XiReactionRate

Description
    Writes the turbulent flame-speed and reaction-rate volScalarFields for the
    Xi-based combustion models.

    Example of function object specification:
    \verbatim
    XiReactionRate
    {
        type        XiReactionRate;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property  | Description                 | Required  | Default value
        type      | type name: XiReactionRate   | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    XiReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_XiReactionRate_H
#define functionObjects_XiReactionRate_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class XiReactionRate Declaration
\*---------------------------------------------------------------------------*/

class XiReactionRate
:
    public fvMeshFunctionObject
{
    // Private member functions

        //- Disallow default bitwise copy construct
        XiReactionRate(const XiReactionRate&);

        //- Disallow default bitwise assignment
        void operator=(const XiReactionRate&);


public:

    //- Runtime type information
    TypeName("XiReactionRate");


    // Constructors

        //- Construct from Time and dictionary
        XiReactionRate
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~XiReactionRate();


    // Member Functions

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
