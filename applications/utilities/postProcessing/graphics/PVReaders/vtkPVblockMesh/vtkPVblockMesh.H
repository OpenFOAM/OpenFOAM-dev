/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::vtkPVblockMesh

Description
    Provides a reader interface for OpenFOAM blockMesh to VTK interaction

SourceFiles
    vtkPVblockMesh.C
    vtkPVblockMeshConvert.C
    vtkPVblockMeshUpdate.C
    vtkPVblockMeshUtils.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVblockMesh_H
#define vtkPVblockMesh_H

// do not include legacy strstream headers
#ifndef  VTK_EXCLUDE_STRSTREAM_HEADERS
    #define VTK_EXCLUDE_STRSTREAM_HEADERS
#endif

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"

#include "primitivePatch.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPVblockMeshReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// OpenFOAM class forward declarations
class argList;
class Time;
class blockMesh;

template<class Type> class List;

/*---------------------------------------------------------------------------*\
                     Class vtkPVblockMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkPVblockMesh
{
    // Private classes

        //- Bookkeeping for GUI checklists and the multi-block organization
        class arrayRange
        {
            const char *name_;
            int block_;
            int start_;
            int size_;

        public:

            arrayRange(const char *name, const int blockNo=0)
            :
                name_(name),
                block_(blockNo),
                start_(0),
                size_(0)
            {}

            //- Return the block holding these datasets
            int block() const
            {
                return block_;
            }

            //- Assign block number, return previous value
            int block(int blockNo)
            {
                int prev = block_;
                block_ = blockNo;
                return prev;
            }

            //- Return block name
            const char* name() const
            {
                return name_;
            }

            //- Return array start index
            int start() const
            {
                return start_;
            }

            //- Return array end index
            int end() const
            {
                return start_ + size_;
            }

            //- Return sublist size
            int size() const
            {
                return size_;
            }

            bool empty() const
            {
                return !size_;
            }

            //- Reset the size to zero and optionally assign a new start
            void reset(const int startAt = 0)
            {
                start_ = startAt;
                size_ = 0;
            }

            //- Increment the size
            void operator+=(const int n)
            {
                size_ += n;
            }
        };


    // Private Data

        //- Access to the controlling vtkPVblockMeshReader
        vtkPVblockMeshReader* reader_;

        //- OpenFOAM time control
        autoPtr<Time> dbPtr_;

        //- OpenFOAM mesh
        blockMesh* meshPtr_;

        //- The mesh region
        word meshRegion_;

        //- The mesh directory for the region
        fileName meshDir_;

        //- Selected geometrical parts
        boolList blockStatus_;

        //- Selected curved edges
        boolList edgeStatus_;

        //- First instance and size of bleckMesh blocks
        //  used to index into blockStatus_
        arrayRange arrayRangeBlocks_;

        //- First instance and size of CurvedEdges (only partially used)
        arrayRange arrayRangeEdges_;

        //- First instance and size of block corners (only partially used)
        arrayRange arrayRangeCorners_;

        //- List of point numbers for rendering to window
        List<vtkTextActor*> pointNumberTextActorsPtrs_;

    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            vtkDataSet* dataset,
            const arrayRange&,
            const label datasetNo,
            const std::string& datasetName
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            const arrayRange&,
            const label datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            const arrayRange&
        );

        //- Update boolList from GUI selection
        static void updateBoolListStatus
        (
            boolList&,
            vtkDataArraySelection*
        );

        //- Reset data counters
        void resetCounters();

        // Update information helper functions

            //- Internal block info
            void updateInfoBlocks(vtkDataArraySelection*);

            //- Block curved edges info
            void updateInfoEdges(vtkDataArraySelection*);

        // Update helper functions

            //- OpenFOAM mesh
            void updateFoamMesh();

        // Mesh conversion functions

            //- Mesh blocks
            void convertMeshBlocks(vtkMultiBlockDataSet*, int& blockNo);

            //- Mesh curved edges
            void convertMeshEdges(vtkMultiBlockDataSet*, int& blockNo);

            //- Mesh corners
            void convertMeshCorners(vtkMultiBlockDataSet*, int& blockNo);


       // GUI selection helper functions

            //- Retrieve the current selections
            static wordHashSet getSelected(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static wordHashSet getSelected
            (
                vtkDataArraySelection*,
                const arrayRange&
            );

            //- Retrieve the current selections
            static stringList getSelectedArrayEntries(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const arrayRange&
            );

            //- Set selection(s)
            static void setSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const stringList&
            );


public:

    //- Static data members

        ClassName("vtkPVblockMesh");


    // Constructors

        //- Construct from components
        vtkPVblockMesh
        (
            const char* const FileName,
            vtkPVblockMeshReader* reader
        );

        //- Disallow default bitwise copy construction
        vtkPVblockMesh(const vtkPVblockMesh&) = delete;


    //- Destructor
    ~vtkPVblockMesh();


    // Member Functions

        //- Update
        void updateInfo();

        void Update(vtkMultiBlockDataSet* output);

        //- Clean any storage
        void CleanUp();

        //- Add/remove point numbers to/from the view
        void renderPointNumbers(vtkRenderer*, const bool show);

     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const vtkPVblockMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
