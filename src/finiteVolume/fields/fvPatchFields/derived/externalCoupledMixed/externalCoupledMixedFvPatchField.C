/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "externalCoupledMixedFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::externalCoupledMixedFvPatchField<Type>::lockName = "OpenFOAM";

template<class Type>
Foam::string
Foam::externalCoupledMixedFvPatchField<Type>::patchKey = "# Patch: ";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::externalCoupledMixedFvPatchField<Type>::baseDir() const
{
    word regionName(this->internalField().mesh().name());
    if (regionName == polyMesh::defaultRegion)
    {
        regionName = ".";
    }

    fileName result(commsDir_/regionName);
    result.clean();

    return result;
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::setMaster
(
    const labelList& patchIDs
)
{
    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->internalField());

    volFieldType& vf = const_cast<volFieldType&>(cvf);

    typename volFieldType::Boundary& bf = vf.boundaryFieldRef();

    // number of patches can be different in parallel...
    label nPatch = bf.size();
    reduce(nPatch, maxOp<label>());

    offsets_.setSize(nPatch);
    forAll(offsets_, i)
    {
        offsets_[i].setSize(Pstream::nProcs());
        offsets_[i] = 0;
    }

    // set the master patch
    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        patchType& pf = refCast<patchType>(bf[patchi]);

        offsets_[patchi][Pstream::myProcNo()] = pf.size();

        if (i == 0)
        {
            pf.master() = true;
        }
        else
        {
            pf.master() = false;
        }
    }

    // set the patch offsets
    int tag = Pstream::msgType() + 1;
    forAll(offsets_, i)
    {
        Pstream::gatherList(offsets_[i], tag);
        Pstream::scatterList(offsets_[i], tag);
    }

    label patchOffset = 0;
    forAll(offsets_, patchi)
    {
        label sumOffset = 0;
        List<label>& procOffsets = offsets_[patchi];

        forAll(procOffsets, proci)
        {
            label o = procOffsets[proci];
            if (o > 0)
            {
                procOffsets[proci] = patchOffset + sumOffset;
                sumOffset += o;
            }
        }
        patchOffset += sumOffset;
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeGeometry
(
    OFstream& osPoints,
    OFstream& osFaces
) const
{
    int tag = Pstream::msgType() + 1;

    const label proci = Pstream::myProcNo();
    const polyPatch& p = this->patch().patch();
    const polyMesh& mesh = p.boundaryMesh().mesh();

    labelList pointToGlobal;
    labelList uniquePointIDs;
    (void)mesh.globalData().mergePoints
    (
        p.meshPoints(),
        p.meshPointMap(),
        pointToGlobal,
        uniquePointIDs
    );

    List<pointField> allPoints(Pstream::nProcs());
    allPoints[proci] = pointField(mesh.points(), uniquePointIDs);
    Pstream::gatherList(allPoints, tag);

    List<faceList> allFaces(Pstream::nProcs());
    faceList& patchFaces = allFaces[proci];
    patchFaces = p.localFaces();
    forAll(patchFaces, facei)
    {
        inplaceRenumber(pointToGlobal, patchFaces[facei]);
    }

    Pstream::gatherList(allFaces, tag);

    if (Pstream::master())
    {
        pointField pts
        (
            ListListOps::combine<pointField>(allPoints, accessOp<pointField>())
        );

        // write points
        osPoints << patchKey.c_str() << this->patch().name() << pts << endl;

        faceList fcs
        (
            ListListOps::combine<faceList>(allFaces, accessOp<faceList>())
        );

        // write faces
        osFaces<< patchKey.c_str() << this->patch().name() << fcs << endl;
    }
}


template<class Type>
Foam::fileName Foam::externalCoupledMixedFvPatchField<Type>::lockFile() const
{
    return fileName(baseDir()/(lockName + ".lock"));
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::createLockFile() const
{
    if (!master_ || !Pstream::master())
    {
        return;
    }

    const fileName fName(lockFile());
    IFstream is(fName);

    // only create lock file if it doesn't already exist
    if (!is.good())
    {
        if (log_)
        {
            Info<< type() << ": creating lock file" << endl;
        }

        OFstream os(fName);
        os  << "lock file";
        os.flush();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::removeLockFile() const
{
    if (!master_ || !Pstream::master())
    {
        return;
    }

    if (log_)
    {
        Info<< type() << ": removing lock file" << endl;
    }

    rm(lockFile());
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::startWait() const
{
    // only wait on master patch

    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->internalField());

    const typename volFieldType::Boundary& bf =
        cvf.boundaryField();

    forAll(coupledPatchIDs_, i)
    {
        label patchi = coupledPatchIDs_[i];

        const patchType& pf = refCast<const patchType>(bf[patchi]);

        if (pf.master())
        {
            pf.wait();
            break;
        }
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::wait() const
{
    const fileName fName(lockFile());
    label found = 0;
    label totalTime = 0;

    if (log_)
    {
        Info<< type() << ": beginning wait for lock file " << fName << endl;
    }

    while (found == 0)
    {
        if (Pstream::master())
        {
            if (totalTime > timeOut_)
            {
                FatalErrorInFunction
                    << "Wait time exceeded time out time of " << timeOut_
                    << " s" << abort(FatalError);
            }

            IFstream is(fName);

            if (is.good())
            {
                found++;

                if (log_)
                {
                    Info<< type() << ": found lock file " << fName << endl;
                }
            }
            else
            {
                sleep(waitInterval_);
                totalTime += waitInterval_;

                if (log_)
                {
                    Info<< type() << ": wait time = " << totalTime << endl;
                }
            }
        }

        // prevent other procs from racing ahead
        reduce(found, sumOp<label>());
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::initialiseRead
(
    IFstream& is
) const
{
    if (!is.good())
    {
        FatalErrorInFunction
            << "Unable to open data transfer file " << is.name()
            << " for patch " << this->patch().name()
            << exit(FatalError);
    }

    label offset = offsets_[this->patch().index()][Pstream::myProcNo()];

    // scan forward to start reading at relevant line/position
    string line;
    for (label i = 0; i < offset; i++)
    {
        if (is.good())
        {
            is.getLine(line);
        }
        else
        {
            FatalErrorInFunction
                << "Unable to scan forward to appropriate read position for "
                << "data transfer file " << is.name()
                << " for patch " << this->patch().name()
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::readData
(
    const fileName& transferFile
)
{
    // read data passed back from external source
    IFstream is(transferFile + ".in");

    // pre-process the input transfer file
    initialiseRead(is);

    // read data from file
    forAll(this->patch(), facei)
    {
        if (is.good())
        {
            is  >> this->refValue()[facei]
                >> this->refGrad()[facei]
                >> this->valueFraction()[facei];
        }
        else
        {
            FatalErrorInFunction
                << "Insufficient data for patch "
                << this->patch().name()
                << " in file " << is.name()
                << exit(FatalError);
        }
    }

    initialised_ = true;

    // update the value from the mixed condition
    mixedFvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeData
(
    const fileName& transferFile
) const
{
    if (!master_)
    {
        return;
    }

    OFstream os(transferFile);

    writeHeader(os);

    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->internalField());

    const typename volFieldType::Boundary& bf =
        cvf.boundaryField();

    forAll(coupledPatchIDs_, i)
    {
        label patchi = coupledPatchIDs_[i];

        const patchType& pf = refCast<const patchType>(bf[patchi]);

        pf.transferData(os);
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeHeader
(
    OFstream& os
) const
{
    os  << "# Values: magSf value snGrad" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    commsDir_("unknown-commsDir"),
    fName_("unknown-fName"),
    waitInterval_(0),
    timeOut_(0),
    calcFrequency_(0),
    initByExternal_(false),
    log_(false),
    master_(false),
    offsets_(),
    initialised_(false),
    coupledPatchIDs_()
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    commsDir_(ptf.commsDir_),
    fName_(ptf.fName_),
    waitInterval_(ptf.waitInterval_),
    timeOut_(ptf.timeOut_),
    calcFrequency_(ptf.calcFrequency_),
    initByExternal_(ptf.initByExternal_),
    log_(ptf.log_),
    master_(ptf.master_),
    offsets_(ptf.offsets_),
    initialised_(ptf.initialised_),
    coupledPatchIDs_(ptf.coupledPatchIDs_)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    commsDir_(dict.lookup("commsDir")),
    fName_(dict.lookup("file")),
    waitInterval_(dict.lookupOrDefault("waitInterval", 1)),
    timeOut_(dict.lookupOrDefault("timeOut", 100*waitInterval_)),
    calcFrequency_(dict.lookupOrDefault("calcFrequency", 1)),
    initByExternal_(readBool(dict.lookup("initByExternal"))),
    log_(dict.lookupOrDefault("log", false)),
    master_(true),
    offsets_(),
    initialised_(false),
    coupledPatchIDs_()
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    commsDir_.expand();

    if (Pstream::master())
    {
        mkDir(baseDir());
    }

    if (!initByExternal_)
    {
        createLockFile();
    }

    // initialise as a fixed value
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf
)
:
    mixedFvPatchField<Type>(ecmpf),
    commsDir_(ecmpf.commsDir_),
    fName_(ecmpf.fName_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    initByExternal_(ecmpf.initByExternal_),
    log_(ecmpf.log_),
    master_(ecmpf.master_),
    offsets_(ecmpf.offsets_),
    initialised_(ecmpf.initialised_),
    coupledPatchIDs_(ecmpf.coupledPatchIDs_)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ecmpf, iF),
    commsDir_(ecmpf.commsDir_),
    fName_(ecmpf.fName_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    initByExternal_(ecmpf.initByExternal_),
    log_(ecmpf.log_),
    master_(ecmpf.master_),
    offsets_(ecmpf.offsets_),
    initialised_(ecmpf.initialised_),
    coupledPatchIDs_(ecmpf.coupledPatchIDs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
~externalCoupledMixedFvPatchField()
{
    removeLockFile();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::initialise
(
    const fileName& transferFile
)
{
    if (initialised_)
    {
        return;
    }

    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->internalField());

    volFieldType& vf = const_cast<volFieldType&>(cvf);

    typename volFieldType::Boundary& bf = vf.boundaryFieldRef();

    // identify all coupled patches
    DynamicList<label> coupledPatchIDs(bf.size());

    forAll(bf, patchi)
    {
        if (isA<patchType>(bf[patchi]))
        {
            coupledPatchIDs.append(patchi);
        }
    }

    coupledPatchIDs_.transfer(coupledPatchIDs);


    // initialise by external solver, or just set the master patch
    if (initByExternal_)
    {
        // remove lock file, signalling external source to execute
//        removeLockFile();

        forAll(coupledPatchIDs_, i)
        {
            label patchi = coupledPatchIDs_[i];

            patchType& pf = refCast<patchType>(bf[patchi]);

            pf.setMaster(coupledPatchIDs_);
        }


        // wait for initial data to be made available
        startWait();

        // read the initial data
        if (master_)
        {
            forAll(coupledPatchIDs_, i)
            {
                label patchi = coupledPatchIDs_[i];

                patchType& pf = refCast<patchType>(bf[patchi]);

                pf.readData(transferFile);
            }
        }
    }
    else
    {
        setMaster(coupledPatchIDs_);
    }

    initialised_ = true;
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes comms
)
{
    if (!initialised_ || this->db().time().timeIndex() % calcFrequency_ == 0)
    {
        const fileName transferFile(baseDir()/fName_);

        // initialise the coupling
        initialise(transferFile);

        // write data for external source
        writeData(transferFile + ".out");

        // remove lock file, signalling external source to execute
        removeLockFile();

        // wait for response
        startWait();

        if (master_ && Pstream::master())
        {
            // remove old data file from OpenFOAM
            rm(transferFile + ".out");
        }

        // read data passed back from external source
        readData(transferFile);

        // create lock file for external source
        createLockFile();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::transferData
(
    OFstream& os
) const
{
    if (log_)
    {
        Info<< type() << ": writing data to " << os.name()  << endl;
    }

    if (Pstream::parRun())
    {
        int tag = Pstream::msgType() + 1;

        List<Field<scalar>> magSfs(Pstream::nProcs());
        magSfs[Pstream::myProcNo()].setSize(this->patch().size());
        magSfs[Pstream::myProcNo()] = this->patch().magSf();
        Pstream::gatherList(magSfs, tag);

        List<Field<Type>> values(Pstream::nProcs());
        values[Pstream::myProcNo()].setSize(this->patch().size());
        values[Pstream::myProcNo()] = this->refValue();
        Pstream::gatherList(values, tag);

        List<Field<Type>> snGrads(Pstream::nProcs());
        snGrads[Pstream::myProcNo()].setSize(this->patch().size());
        snGrads[Pstream::myProcNo()] = this->snGrad();
        Pstream::gatherList(snGrads, tag);

        if (Pstream::master())
        {
            forAll(values, proci)
            {
                const Field<scalar>& magSf = magSfs[proci];
                const Field<Type>& value = values[proci];
                const Field<Type>& snGrad = snGrads[proci];

                forAll(magSf, facei)
                {
                    os  << magSf[facei] << token::SPACE
                        << value[facei] << token::SPACE
                        << snGrad[facei] << nl;
                }
            }

            os.flush();
        }
    }
    else
    {
        const Field<scalar>& magSf(this->patch().magSf());
        const Field<Type>& value(this->refValue());
        const Field<Type> snGrad(this->snGrad());

        forAll(magSf, facei)
        {
            os  << magSf[facei] << token::SPACE
                << value[facei] << token::SPACE
                << snGrad[facei] << nl;
        }

        os.flush();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeGeometry() const
{
    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->internalField());

    const typename volFieldType::Boundary& bf =
        cvf.boundaryField();

    OFstream osPoints(baseDir()/"patchPoints");
    OFstream osFaces(baseDir()/"patchFaces");

    if (log_)
    {
        Info<< "writing collated patch points to: " << osPoints.name() << nl
            << "writing collated patch faces to: " << osFaces.name() << endl;
    }

    forAll(bf, patchi)
    {
        if (isA<patchType>(bf[patchi]))
        {
            const patchType& pf = refCast<const patchType>(bf[patchi]);

            pf.writeGeometry(osPoints, osFaces);
        }
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::write(Ostream& os) const
{
    mixedFvPatchField<Type>::write(os);

    writeEntry(os, "commsDir", commsDir_);
    writeEntry(os, "file", fName_);
    writeEntry(os, "waitInterval", waitInterval_);
    writeEntry(os, "timeOut", timeOut_);
    writeEntry(os, "calcFrequency", calcFrequency_);
    writeEntry(os, "initByExternal", initByExternal_);
    writeEntry(os, "log", log_);

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
