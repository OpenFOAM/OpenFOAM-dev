/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::functionObjects::wallHeatFlux

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'wallHeatFlux'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    wallHeatFlux1
    {
        type        wallHeatFlux;
        libs        ("libfieldFunctionObjects.so");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFlux    | yes        |
        patches  | list of patches to process | no         | all wall patches
        region   | region to be evaluated     | no         | default region
    \endtable

    Note:
        Writing field 'wallHeatFlux' is done by default, but it can be
        overridden by defining an empty \c objects list. For details see
        writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFlux_H
#define functionObjects_wallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "HashSet.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class wallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFlux
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Data

        //- The name of the phase
        word phaseName_;

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the heat-flux
        tmp<volScalarField> calcWallHeatFlux(const surfaceScalarField& q);


public:

    //- Runtime type information
    TypeName("wallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        wallHeatFlux(const wallHeatFlux&) = delete;


    //- Destructor
    virtual ~wallHeatFlux();


    // Member Functions

        //- Read the wallHeatFlux data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallHeatFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
