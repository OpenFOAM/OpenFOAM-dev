/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::vtkPV4blockMesh

Description
    Provides a reader interface for OpenFOAM blockMesh to VTK interaction

SourceFiles
    vtkPV4blockMesh.C
    vtkPV4blockMeshConvert.C
    vtkPV4blockMeshUpdate.C
    vtkPV4blockMeshUtils.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV4blockMesh_H
#define vtkPV4blockMesh_H

// do not include legacy strstream headers
#ifndef  VTK_EXCLUDE_STRSTREAM_HEADERS
# define VTK_EXCLUDE_STRSTREAM_HEADERS
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
class vtkPV4blockMeshReader;
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
                     Class vtkPV4blockMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkPV4blockMesh
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

        //- Access to the controlling vtkPV4blockMeshReader
        vtkPV4blockMeshReader* reader_;

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

            //- block curved edges info
            void updateInfoEdges(vtkDataArraySelection*);

        // Update helper functions

            //- OpenFOAM mesh
            void updateFoamMesh();

        // Mesh conversion functions

            //- mesh blocks
            void convertMeshBlocks(vtkMultiBlockDataSet*, int& blockNo);

            //- mesh curved edges
            void convertMeshEdges(vtkMultiBlockDataSet*, int& blockNo);

            //- mesh corners
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


        //- Disallow default bitwise copy construct
        vtkPV4blockMesh(const vtkPV4blockMesh&);

        //- Disallow default bitwise assignment
        void operator=(const vtkPV4blockMesh&);


public:

    //- Static data members

        ClassName("vtkPV4blockMesh");


    // Constructors

        //- Construct from components
        vtkPV4blockMesh
        (
            const char* const FileName,
            vtkPV4blockMeshReader* reader
        );


    //- Destructor
    ~vtkPV4blockMesh();


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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
