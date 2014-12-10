/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::vtkPV3Readers

Description
    A collection of helper functions when building a reader interface in
    ParaView3.

SourceFiles
    vtkPV3Readers.C

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3Readers_H
#define vtkPV3Readers_H

// do not include legacy strstream headers
#ifndef  VTK_EXCLUDE_STRSTREAM_HEADERS
# define VTK_EXCLUDE_STRSTREAM_HEADERS
#endif

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"
#include "HashSet.H"


// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPV3FoamReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtkPV3Readers
{
    //- Declare name of the class and its debug switch
    NamespaceName("vtkPV3Readers");

    //- Bookkeeping for GUI checklists and the multi-block organization
    class partInfo
    {
        const char *name_;
        int block_;
        int start_;
        int size_;

    public:

        partInfo(const char *name, const int blockNo=0)
        :
            name_(name),
            block_(blockNo),
            start_(-1),
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

        const char* name() const
        {
            return name_;
        }

        int start() const
        {
            return start_;
        }

        int end() const
        {
            return start_ + size_;
        }

        int size() const
        {
            return size_;
        }

        bool empty() const
        {
            return !size_;
        }

        void reset()
        {
            start_ = -1;
            size_ = 0;
        }

        //- Assign new start and reset the size
        void operator=(const int i)
        {
            start_ = i;
            size_ = 0;
        }

        //- Increment the size
        void operator+=(const int n)
        {
            size_ += n;
        }
    };


    //- Convenience method use to convert the readers from VTK 5
    //  multiblock API to the current composite data infrastructure
    void AddToBlock
    (
        vtkMultiBlockDataSet* output,
        vtkDataSet* dataset,
        const partInfo& selector,
        const label datasetNo,
        const std::string& datasetName
    );


    //- Convenience method use to convert the readers from VTK 5
    // multiblock API to the current composite data infrastructure
    vtkDataSet* GetDataSetFromBlock
    (
        vtkMultiBlockDataSet* output,
        const partInfo& selector,
        const label datasetNo
    );

    //- Convenience method use to convert the readers from VTK 5
    // multiblock API to the current composite data infrastructure
    // ununsed at the moment
    label GetNumberOfDataSets
    (
        vtkMultiBlockDataSet* output,
        const partInfo& selector
    );


    //- Retrieve the current selections as a wordHashSet
    wordHashSet getSelected
    (
        vtkDataArraySelection* select
    );


    //- Retrieve a sub-list of the current selections
    wordHashSet getSelected
    (
        vtkDataArraySelection*,
        const partInfo&
    );


    //- Retrieve the current selections
    stringList getSelectedArrayEntries(vtkDataArraySelection*);

    //- Retrieve a sub-list of the current selections
    stringList getSelectedArrayEntries
    (
        vtkDataArraySelection* select,
        const partInfo& selector
    );


    //- Set selection(s)
    void setSelectedArrayEntries
    (
        vtkDataArraySelection*,
        const stringList&
    );



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtkPV3

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
