/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation.H"
#include "regionCoupledBaseGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionCoupledBaseGAMGInterface::regionCoupledBaseGAMGInterface
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
    fineRegionCoupledLduInterface_
    (
        refCast<const regionCoupledLduInterface>(fineInterface)
    )
{
    // Construct face agglomeration from cell agglomeration
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        forAll(localRestrictAddressing, ffi)
        {
            label curMaster = localRestrictAddressing[ffi];

            Map<label>::const_iterator fnd = masterToCoarseFace.find
            (
                curMaster
            );

            if (fnd == masterToCoarseFace.end())
            {
                // New coarse face
                label coarseI = dynFaceCells.size();
                dynFaceRestrictAddressing.append(coarseI);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseI);
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

    /*
    // On the owner side construct the AMI
    if (fineRegionCoupledLduInterface_.owner())
    {
        const polyMesh& nbrMesh =
            fineRegionCoupledLduInterface_.nbrMesh();

        if
        (
            nbrMesh.foundObject<GAMGAgglomeration>(GAMGAgglomeration::typeName)
        )
        {
            const GAMGAgglomeration& nbrAgg = nbrMesh.thisDb().lookupObject
            <
                GAMGAgglomeration
            >
            (
                GAMGAgglomeration::typeName
            );

            label nbrLevel(-1);
            if (nbrAgg.size() > fineLevelIndex)
            {
                nbrLevel = fineLevelIndex;
            }
            else
            {
                nbrLevel = nbrAgg.size() - 1;
            }

            const labelField& nbrRestrictMap =
                nbrAgg.restrictAddressing(nbrLevel);

            const labelUList& nbrFaceCells =
                nbrLduInterface
                (
                    nbrLevel,
                    neighbPatchID()
                ).faceCells();


            const IndirectList<label> nbrPatchRestrictMap
            (
                nbrRestrictMap,
                nbrFaceCells
            );

            labelList nbrFaceRestrictAddressing;
            {
                // From face to coarse face
                DynamicList<label> dynNbrFaceRestrictAddressing
                (
                    nbrPatchRestrictMap.size()
                );

                Map<label> masterToCoarseFace(nbrPatchRestrictMap.size());

                forAll(nbrPatchRestrictMap, ffi)
                {
                    label curMaster = nbrPatchRestrictMap[ffi];

                    Map<label>::const_iterator fnd = masterToCoarseFace.find
                    (
                        curMaster
                    );

                    if (fnd == masterToCoarseFace.end())
                    {
                        // New coarse face
                        label coarseI = masterToCoarseFace.size();
                        dynNbrFaceRestrictAddressing.append(coarseI);
                        masterToCoarseFace.insert(curMaster, coarseI);
                    }
                    else
                    {
                        // Already have coarse face
                        dynNbrFaceRestrictAddressing.append(fnd());
                    }
                }

                nbrFaceRestrictAddressing.transfer
                (
                    dynNbrFaceRestrictAddressing
                );
            }

            amiPtr_.reset
            (
                new AMIPatchToPatchInterpolation
                (
                    fineRegionCoupledLduInterface_.AMI(),
                    faceRestrictAddressing_,
                    nbrFaceRestrictAddressing
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << " GAMGAgglomeration was not found in the nbr mesh. "
                << " Check on the cacheAgglomeration flag in fvSolution"
                << exit(FatalError);
        }
    }
    */

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledBaseGAMGInterface::~regionCoupledBaseGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::regionCoupledBaseGAMGInterface::
internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    /*
    //const labelUList& nbrFaceCells = neighbPatch().faceCells();

    const labelUList& nbrFaceCells = nbrLduInterface().faceCells();

    tmp<labelField> tpnf(new labelField(nbrFaceCells.size()));
    labelField& pnf = tpnf();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }
    */
    tmp<labelField> tpnf(new labelField(iF));

    return tpnf;
}


// ************************************************************************* //
