/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "processorGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "HashTable.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        processorGAMGInterface,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        processorGAMGInterface,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorGAMGInterface::processorGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces
    ),
    comm_(coarseComm),
    myProcNo_(refCast<const processorLduInterface>(fineInterface).myProcNo()),
    neighbProcNo_
    (
        refCast<const processorLduInterface>(fineInterface).neighbProcNo()
    ),
    transform_
    (
        refCast<const processorLduInterface>(fineInterface).transform()
    ),
    tag_(refCast<const processorLduInterface>(fineInterface).tag())
{
    // From coarse face to coarse cell
    DynamicList<label> dynFaceCells(localRestrictAddressing.size());
    // From fine face to coarse face
    DynamicList<label> dynFaceRestrictAddressing
    (
        localRestrictAddressing.size()
    );

    // From coarse cell pair to coarse face
    HashTable<label, labelPair, labelPair::Hash<>> cellsToCoarseFace
    (
        2*localRestrictAddressing.size()
    );

    forAll(localRestrictAddressing, ffi)
    {
        labelPair cellPair;

        // Do switching on master/slave indexes based on the owner/neighbour of
        // the processor index such that both sides get the same answer.
        if (myProcNo() < neighbProcNo())
        {
            // Master side
            cellPair = labelPair
            (
                localRestrictAddressing[ffi],
                neighbourRestrictAddressing[ffi]
            );
        }
        else
        {
            // Slave side
            cellPair = labelPair
            (
                neighbourRestrictAddressing[ffi],
                localRestrictAddressing[ffi]
            );
        }

        HashTable<label, labelPair, labelPair::Hash<>>::const_iterator fnd =
            cellsToCoarseFace.find(cellPair);

        if (fnd == cellsToCoarseFace.end())
        {
            // New coarse face
            label coarseI = dynFaceCells.size();
            dynFaceRestrictAddressing.append(coarseI);
            dynFaceCells.append(localRestrictAddressing[ffi]);
            cellsToCoarseFace.insert(cellPair, coarseI);
        }
        else
        {
            // Already have coarse face
            dynFaceRestrictAddressing.append(fnd());
        }
    }

    faceCells_.transfer(dynFaceCells);
    faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
}


Foam::processorGAMGInterface::processorGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const labelUList& faceCells,
    const labelUList& faceRestrictAddressing,
    const label coarseComm,
    const label myProcNo,
    const label neighbProcNo,
    const transformer& transform,
    const int tag
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces,
        faceCells,
        faceRestrictAddressing
    ),
    comm_(coarseComm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    transform_(transform),
    tag_(tag)
{}


Foam::processorGAMGInterface::processorGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    GAMGInterface(index, coarseInterfaces, is),
    comm_(readLabel(is)),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is)),
    transform_(is),
    tag_(readLabel(is))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::processorGAMGInterface::~processorGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, interfaceInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    tmp<Field<label>> tfld(receive<label>(commsType, this->size()));

    return tfld;
}


void Foam::processorGAMGInterface::write(Ostream& os) const
{
    GAMGInterface::write(os);
    os  << token::SPACE << comm_
        << token::SPACE << myProcNo_
        << token::SPACE << neighbProcNo_
        << token::SPACE << transform_
        << token::SPACE << tag_;
}


// ************************************************************************* //
