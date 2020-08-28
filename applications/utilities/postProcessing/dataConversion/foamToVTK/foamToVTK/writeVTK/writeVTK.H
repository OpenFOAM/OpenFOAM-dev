/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::functionObjects::writeVTK

Description
    This functionObject writes objects registered to the database in VTK format
    using the foamToVTK library.

    Currently only the writing of the cell-values of volFields is supported but
    support for other field types, patch fields, Lagrangian data etc. will be
    added.

    Example of function object specification:
    \verbatim
        writeVTK1
        {
            type        writeVTK;
            libs ("libutilityFunctionObjects.so");
            ...
            objectNames (obj1 obj2);
        }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: writeVTK     | yes         |
        objectNames  | objects to write        | yes         |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    writeVTK.C
    IOwriteVTK.H

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeVTK_H
#define functionObjects_writeVTK_H

#include "fvMeshFunctionObject.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class writeVTK Declaration
\*---------------------------------------------------------------------------*/

class writeVTK
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Names of objects
        wordReList objectNames_;


    // Private Member Functions

        template<class GeoField>
        UPtrList<const GeoField> lookupFields() const;


public:

    //- Runtime type information
    TypeName("writeVTK");


    // Constructors

        //- Construct from Time and dictionary
        writeVTK
        (
            const word& name,
            const Time& t,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        writeVTK(const writeVTK&) = delete;


    //- Destructor
    virtual ~writeVTK();


    // Member Functions

        //- Read the writeVTK data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the writeVTK
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const writeVTK&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "writeVTKTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
