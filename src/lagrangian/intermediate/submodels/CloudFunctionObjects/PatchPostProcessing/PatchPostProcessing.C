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

\*---------------------------------------------------------------------------*/

#include "PatchPostProcessing.H"
#include "Pstream.H"
#include "stringListOps.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchPostProcessing<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchi)
        {
            return i;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::write()
{
    forAll(patchData_, i)
    {
        List<List<scalar>> procTimes(Pstream::nProcs());
        procTimes[Pstream::myProcNo()] = times_[i];
        Pstream::gatherList(procTimes);

        List<List<string>> procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = patchData_[i];
        Pstream::gatherList(procData);

        if (Pstream::master())
        {
            const fvMesh& mesh = this->owner().mesh();

            // Create directory if it doesn't exist
            mkDir(this->writeTimeDir());

            const word& patchName = mesh.boundaryMesh()[patchIDs_[i]].name();

            OFstream patchOutFile
            (
                this->writeTimeDir()/patchName + ".post",
                IOstream::ASCII,
                IOstream::currentVersion,
                mesh.time().writeCompression()
            );

            List<string> globalData;
            globalData = ListListOps::combine<List<string>>
            (
                procData,
                accessOp<List<string>>()
            );

            List<scalar> globalTimes;
            globalTimes = ListListOps::combine<List<scalar>>
            (
                procTimes,
                accessOp<List<scalar>>()
            );

            labelList indices;
            sortedOrder(globalTimes, indices);

            string header("# Time currentProc " + parcelType::propertyList_);
            patchOutFile<< header.c_str() << nl;

            forAll(globalTimes, i)
            {
                label dataI = indices[i];

                patchOutFile
                    << globalTimes[dataI] << ' '
                    << globalData[dataI].c_str()
                    << nl;
            }
        }

        patchData_[i].clearStorage();
        times_[i].clearStorage();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    maxStoredParcels_(readScalar(this->coeffDict().lookup("maxStoredParcels"))),
    patchIDs_(),
    times_(),
    patchData_()
{
    const wordList allPatchNames = owner.mesh().boundaryMesh().names();
    wordList patchName(this->coeffDict().lookup("patches"));

    labelHashSet uniquePatchIDs;
    forAllReverse(patchName, i)
    {
        labelList patchIDs = findStrings(patchName[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching " << patchName[i]
                << endl;
        }

        uniquePatchIDs.insert(patchIDs);
    }

    patchIDs_ = uniquePatchIDs.toc();

    if (debug)
    {
        forAll(patchIDs_, i)
        {
            const label patchi = patchIDs_[i];
            const word& patchName = owner.mesh().boundaryMesh()[patchi].name();
            Info<< "Post-process patch " << patchName << endl;
        }
    }

    patchData_.setSize(patchIDs_.size());
    times_.setSize(patchIDs_.size());
}


template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const PatchPostProcessing<CloudType>& ppm
)
:
    CloudFunctionObject<CloudType>(ppm),
    maxStoredParcels_(ppm.maxStoredParcels_),
    patchIDs_(ppm.patchIDs_),
    times_(ppm.times_),
    patchData_(ppm.patchData_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::~PatchPostProcessing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();
    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1 && patchData_[localPatchi].size() < maxStoredParcels_)
    {
        times_[localPatchi].append(this->owner().time().value());

        OStringStream data;
        data<< Pstream::myProcNo() << ' ' << p;

        patchData_[localPatchi].append(data.str());
    }
}


// ************************************************************************* //
