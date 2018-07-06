/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::functionObjects::interfaceHeight

Description
    This function object reports the height of the interface above a set of
    locations. For each location, it writes the vertical distance of the
    interface above both the location and the lowest boundary. It also writes
    the point on the interface from which these heights are computed. It uses
    an integral approach, so if there are multiple interfaces above or below a
    location then this method will generate average values.

    Example of function object specification:
    \verbatim
    interfaceHeight1
    {
        type           interfaceHeight;
        libs           ("libfieldFunctionObjects.so");
        alpha          alpha.water;
        locations      ((0 0 0) (10 0 0) (20 0 0));
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name                 | yes      |
        alpha        | name of the alpha field   | no       | alpha
        locations    | list of locations to report the height at | yes |
        liquid       | is the alpha field that of the liquid | no | true
    \endtable

SourceFiles
    interfaceHeight.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceHeight_H
#define interfaceHeight_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class interfaceHeight Declaration
\*---------------------------------------------------------------------------*/

class interfaceHeight
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private data

        //- Name of the alpha field
        word alphaName_;

        //- Is the alpha field that of the liquid under the wave?
        bool liquid_;

        //- List of locations to report the height at
        List<point> locations_;

        //- Interpolation scheme
        word interpolationScheme_;


    // Private Member Functions

        //- Output positions
        void writePositions();


    // Private Enumerations

        //- File enumeration
        enum fileID
        {
            HEIGHT_FILE = 0,
            POSITION_FILE = 1
        };


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("interfaceHeight");


    // Constructors

        //- Construct from Time and dictionary
        interfaceHeight
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interfaceHeight();


    // Member Functions

        //- Read
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
