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

Description
    Misc helper methods and utilities

\*---------------------------------------------------------------------------*/

#include "vtkPV4Readers.H"

// OpenFOAM includes
#include "IFstream.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkInformation.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(vtkPV4Readers, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //! \cond fileScope
    //  Extract up to the first non-word characters
    inline word getFirstWord(const char* str)
    {
        if (str)
        {
            label n = 0;
            while (str[n] && word::valid(str[n]))
            {
                ++n;
            }
            return word(str, n, true);
        }
        else
        {
            return word::null;
        }

    }
    //! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkPV4Readers::AddToBlock
(
    vtkMultiBlockDataSet* output,
    vtkDataSet* dataset,
    const partInfo& selector,
    const label datasetNo,
    const std::string& datasetName
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);

    if (!block)
    {
        if (blockDO)
        {
            FatalErrorIn("Foam::vtkPV4Readers::AddToBlock")
                << "Block already has a vtkDataSet assigned to it"
                << endl;
            return;
        }

        block = vtkMultiBlockDataSet::New();
        output->SetBlock(blockNo, block);
        block->Delete();
    }

    if (debug)
    {
        Info<< "block[" << blockNo << "] has "
            << block->GetNumberOfBlocks()
            <<  " datasets prior to adding set " << datasetNo
            <<  " with name: " << datasetName << endl;
    }

    block->SetBlock(datasetNo, dataset);

    // name the block when assigning dataset 0
    if (datasetNo == 0)
    {
        output->GetMetaData(blockNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            selector.name()
        );
    }

    if (datasetName.size())
    {
        block->GetMetaData(datasetNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            datasetName.c_str()
        );
    }
}


vtkDataSet* Foam::vtkPV4Readers::GetDataSetFromBlock
(
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);

    if (block)
    {
        return vtkDataSet::SafeDownCast(block->GetBlock(datasetNo));
    }

    return 0;
}


// ununsed at the moment
Foam::label Foam::vtkPV4Readers::GetNumberOfDataSets
(
    vtkMultiBlockDataSet* output,
    const partInfo& selector
)
{
    const int blockNo = selector.block();

    vtkDataObject* blockDO = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(blockDO);
    if (block)
    {
        return block->GetNumberOfBlocks();
    }

    return 0;
}


// Foam::word Foam::vtkPV4Readers::getPartName(int partId)
// {
//     return getFirstWord(reader_->GetPartArrayName(partId));
// }


Foam::wordHashSet Foam::vtkPV4Readers::getSelected
(
    vtkDataArraySelection* select
)
{
    int nElem = select->GetNumberOfArrays();
    wordHashSet selections(2*nElem);

    for (int elemI=0; elemI < nElem; ++elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            selections.insert(getFirstWord(select->GetArrayName(elemI)));
        }
    }

    return selections;
}


Foam::wordHashSet Foam::vtkPV4Readers::getSelected
(
    vtkDataArraySelection* select,
    const partInfo& selector
)
{
    int nElem = select->GetNumberOfArrays();
    wordHashSet selections(2*nElem);

    for (int elemI = selector.start(); elemI < selector.end(); ++elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            selections.insert(getFirstWord(select->GetArrayName(elemI)));
        }
    }

    return selections;
}


Foam::stringList Foam::vtkPV4Readers::getSelectedArrayEntries
(
    vtkDataArraySelection* select
)
{
    stringList selections(select->GetNumberOfArrays());
    label nElem = 0;

    forAll(selections, elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            selections[nElem++] = select->GetArrayName(elemI);
        }
    }
    selections.setSize(nElem);


    if (debug)
    {
        label nElem = select->GetNumberOfArrays();
        Info<< "available(";
        for (int elemI = 0; elemI < nElem; ++elemI)
        {
            Info<< " \"" << select->GetArrayName(elemI) << "\"";
        }
        Info<< " )\nselected(";

        forAll(selections, elemI)
        {
            Info<< " " << selections[elemI];
        }
        Info<< " )\n";
    }

    return selections;
}


Foam::stringList Foam::vtkPV4Readers::getSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const partInfo& selector
)
{
    stringList selections(selector.size());
    label nElem = 0;

    for (int elemI = selector.start(); elemI < selector.end(); ++elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            selections[nElem++] = select->GetArrayName(elemI);
        }
    }
    selections.setSize(nElem);


    if (debug)
    {
        Info<< "available(";
        for (int elemI = selector.start(); elemI < selector.end(); ++elemI)
        {
            Info<< " \"" << select->GetArrayName(elemI) << "\"";
        }
        Info<< " )\nselected(";

        forAll(selections, elemI)
        {
            Info<< " " << selections[elemI];
        }
        Info<< " )\n";
    }

    return selections;
}


void Foam::vtkPV4Readers::setSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const stringList& selections
)
{
    const int nElem = select->GetNumberOfArrays();
    select->DisableAllArrays();

    // Loop through entries, setting values from selectedEntries
    for (int elemI=0; elemI < nElem; ++elemI)
    {
        string arrayName(select->GetArrayName(elemI));

        forAll(selections, elemI)
        {
            if (selections[elemI] == arrayName)
            {
                select->EnableArray(arrayName.c_str());
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
