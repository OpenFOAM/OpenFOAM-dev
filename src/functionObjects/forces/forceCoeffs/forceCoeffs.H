/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::functionObjects::forceCoeffs

Description
    Extends the forces functionObject by providing lift, drag and moment
    coefficients.  The data can optionally be output into bins, defined in a
    given direction.

    Example of function object specification:
    \verbatim
    forceCoeffs1
    {
        type        forceCoeffs;
        libs        ("libforces.so");
        ...
        log         yes;
        patches     (walls);
        liftDir     (0 1 0);
        dragDir     (-1 0 0);
        pitchAxis   (0 0 1);
        magUInf     100;
        lRef        3.5;
        Aref        2.2;

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: forceCoeffs  | yes         |
        log          | write force data to standard output | no | no
        patches      | patches included in the forces calculation | yes |
        liftDir      | lift direction          | yes         |
        dragDir      | drag direction          | yes         |
        pitchAxis    | picth axis              | yes         |
        magUInf      | free stream velocity magnitude | yes  |
        lRef         | reference length scale for moment calculations | yes |
        Aref         | reference area          | yes |
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces

SourceFiles
    forceCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceCoeffs_H
#define functionObjects_forceCoeffs_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffs Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffs
:
    public forces
{
    // Private data

        // Force coefficient geometry

            //- Lift
            vector liftDir_;

            //- Drag
            vector  dragDir_;

            //- Pitch
            vector pitchAxis_;


        // Free-stream conditions

            //- Velocity magnitude
            scalar magUInf_;


        // Reference scales

            //- Length
            scalar lRef_;

            //- Area
            scalar Aref_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        forceCoeffs(const forceCoeffs&);

        //- Disallow default bitwise assignment
        void operator=(const forceCoeffs&);


protected:

    //- Output file header information
    virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("forceCoeffs");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffs
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~forceCoeffs();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
