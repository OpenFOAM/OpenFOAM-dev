/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::bXiProgress

Description
    Writes the combustion progress of the Weller b-Xi combustion models.

    Example of function object specification:
    \verbatim
    bXiProgress
    {
        type        bXiProgress;
        libs        ("libfieldFunctionObjects.so");
    }
    \endverbatim

Usage
    \table
        Property  | Description                 | Required  | Default value
        type      | type name: bXiProgress      | yes       |
    \endtable

    Or, using the standard configuration:
    \verbatim
    #includeFunc bXiProgress
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    bXiProgress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_bXiProgress_H
#define functionObjects_bXiProgress_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class bXiProgress Declaration
\*---------------------------------------------------------------------------*/

class bXiProgress
:
    public fvMeshFunctionObject,
    public logFiles
{
protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("bXiProgress");


    // Constructors

        //- Construct from Time and dictionary
        bXiProgress
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bXiProgress(const bXiProgress&) = delete;


    //- Destructor
    virtual ~bXiProgress();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do nothing
        virtual bool execute();

        //- Write the cell-centre fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bXiProgress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
