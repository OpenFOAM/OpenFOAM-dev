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
    Foam::functionObjects::histogram

Description
    Write the volume-weighted histogram of a volScalarField.

    Example:
    \verbatim
    histogram1
    {
        type            histogram;

        libs            ("libfieldFunctionObjects.so");

        field           p;
        cellZone        all;
        nBins           100;
        min             -5;
        max             5;
        setFormat       gnuplot;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: histogram    | yes         |
        field        | Field to analyse        | yes         |
        cellZone     | cellZone                | yes         |
        nBins        | Number of bins for the histogram | yes|
        max          | Maximum value sampled   | yes         |
        min          | minimum value sampled   | no          | 0
        setFormat    | Output format           | yes         |
    \endtable

    Where \c cellZone options are:
    \plaintable
        cellZone \<name\>  | Looks-up the named cellZone
        cellZone {type \<zoneGeneratorType\>;...} | \\
            Generates the cellZone locally
        cellZone all       | Selects all cells
    \endplaintable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    histogram.C

\*---------------------------------------------------------------------------*/

#ifndef histogram_functionObject_H
#define histogram_functionObject_H

#include "fvMeshFunctionObject.H"
#include "generatedCellZone.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class setWriter;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class histogram Declaration
\*---------------------------------------------------------------------------*/

class histogram
:
    public fvMeshFunctionObject
{
    // Private Data

        //- The cellZone to operate on
        generatedCellZone zone_;

        writeFile file_;

        //- Name of field
        word fieldName_;

        //- Maximum value
        scalar max_;

        //- Minimum value
        scalar min_;

        //- Number of bins
        label nBins_;

        //- Output formatter to write
        autoPtr<setWriter> formatterPtr_;


public:

    //- Runtime type information
    TypeName("histogram");


    // Constructors

        //- Construct from Time and dictionary
        histogram
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        histogram(const histogram&) = delete;


    // Destructor
    virtual ~histogram();


    // Member Functions

        //- Read the histogram data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the histogram and write.
        //  foamPostProcess overrides the usual writeControl behaviour and
        //  forces writing always (used in post-processing mode)
        virtual bool write();


        // Mesh changes

            //- Update for mesh motion
            virtual void movePoints(const polyMesh&);

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const histogram&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
