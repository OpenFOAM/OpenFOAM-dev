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

Namespace
    Foam::meshWriters

Description
    A namespace for holding various types of mesh writers.


Class
    Foam::meshWriter

Description
    write OpenFOAM meshes and/or results to another CFD format
    - currently just STAR-CD

\par Files

    "constant/boundaryRegion" is an IOMap<dictionary> that contains
    the boundary type and names. eg,
    \verbatim
        (
            0
            {
                BoundaryType    wall;
                Label           Default_Boundary_Region;
            }

            1
            {
                BoundaryType    inlet;
                Label           inlet_1;
            }

            ...

            4
            {
                BoundaryType    pressure;
                Label           outlet;
            }
        )
    \endverbatim


SourceFiles
    meshWriterI.H
    meshWriter.C
    meshWriterIO.C

\*---------------------------------------------------------------------------*/

#ifndef meshWriter_H
#define meshWriter_H

#include "polyMesh.H"
#include "boundaryRegion.H"
#include "cellTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshWriter Declaration
\*---------------------------------------------------------------------------*/

class meshWriter
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshWriter(const meshWriter&);

        //- Disallow default bitwise assignment
        void operator=(const meshWriter&);


protected:

    // Protected data

        //- Mesh reference
        const polyMesh& mesh_;

        //- Scaling factor for points (eg, [m] -> [mm])
        scalar scaleFactor_;

        //- Write bnd file
        bool writeBoundary_;

        //- boundaryRegion persistent data saved as a dictionary
        boundaryRegion boundaryRegion_;

        //- cellTable persistent data saved as a dictionary
        cellTable cellTable_;

        //- cellTable IDs for each cell
        labelList cellTableId_;

        //- Pointers to cell shape models
        static const cellModel* unknownModel;
        static const cellModel* tetModel;
        static const cellModel* pyrModel;
        static const cellModel* prismModel;
        static const cellModel* hexModel;


public:

    // Static data members

        //- Specify a default mesh name
        static string defaultMeshName;

    // Constructors

        //- Create a writer obejct
        meshWriter
        (
            const polyMesh&,
            const scalar scaleFactor = 1.0
        );


    //- Destructor
    virtual ~meshWriter();


    // Member Functions

        // Edit

            //- Set points scaling
            void scaleFactor(const scalar scaling)
            {
                scaleFactor_ = scaling;
            }

            //- Suppress writing bnd file
            void noBoundary()
            {
                writeBoundary_ = false;
            }

        // Write

            //- Write volume mesh. Subclass must supply this method
            virtual bool write
            (
                const fileName& timeName = fileName::null
            ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
