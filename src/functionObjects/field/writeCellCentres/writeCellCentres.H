/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
    Foam::functionObjects::writeCellCentres

Description
    Writes the cell-centres volVectorField and the three component fields as
    volScalarFields.

    Example of function object specification:
    \verbatim
    writeCellCentres
    {
        type        writeCellCentres;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property  | Description                 | Required    | Default value
        type      | type name: writeCellCentres | yes       |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    writeCellCentres.C

\*---------------------------------------------------------------------------*/

#ifndef writeCellCentres_functionObject_H
#define writeCellCentres_functionObject_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class writeCellCentres Declaration
\*---------------------------------------------------------------------------*/

class writeCellCentres
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("writeCellCentres");


    // Constructors

        //- Construct from Time and dictionary
        writeCellCentres
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        writeCellCentres(const writeCellCentres&) = delete;


    //- Destructor
    virtual ~writeCellCentres();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeCellCentres&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
