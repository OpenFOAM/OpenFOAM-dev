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

#include "Cloud.H"
#include "Time.H"
#include "IOPosition.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::word Foam::Cloud<ParticleType>::cloudPropertiesName("cloudProperties");


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudUniformProperties()
{
    IOobject dictObj
    (
        cloudPropertiesName,
        time().timeName(),
        "uniform"/cloud::prefix/name(),
        db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    if (dictObj.typeHeaderOk<IOdictionary>(true))
    {
        const IOdictionary uniformPropsDict(dictObj);

        const word procName("processor" + Foam::name(Pstream::myProcNo()));
        if (uniformPropsDict.found(procName))
        {
            uniformPropsDict.subDict(procName).lookup("particleCount")
                >> ParticleType::particleCount_;
        }
    }
    else
    {
        ParticleType::particleCount_ = 0;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeCloudUniformProperties() const
{
    IOdictionary uniformPropsDict
    (
        IOobject
        (
            cloudPropertiesName,
            time().timeName(),
            "uniform"/cloud::prefix/name(),
            db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList np(Pstream::nProcs(), 0);
    np[Pstream::myProcNo()] = ParticleType::particleCount_;

    Pstream::listCombineGather(np, maxEqOp<label>());
    Pstream::listCombineScatter(np);

    forAll(np, i)
    {
        word procName("processor" + Foam::name(i));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add("particleCount", np[i]);
    }

    uniformPropsDict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        time().writeCompression(),
        true
    );
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloud(const bool checkClass)
{
    readCloudUniformProperties();

    IOPosition<Cloud<ParticleType>> ioP(*this);

    bool valid = ioP.headerOk();
    Istream& is = ioP.readStream(checkClass ? typeName : "", valid);
    if (valid)
    {
        ioP.readData(is, *this);
        ioP.close();
    }

    if (!valid && debug)
    {
        Pout<< "Cannot read particle positions file:" << nl
            << "    " << ioP.objectPath() << nl
            << "Assuming the initial cloud contains 0 particles." << endl;
    }

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const bool checkClass
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    globalPositionsPtr_()
{
    checkPatches();

    initCloud(checkClass);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOobject Foam::Cloud<ParticleType>::fieldIOobject
(
    const word& fieldName,
    const IOobject::readOption r
) const
{
    return IOobject
    (
        fieldName,
        time().timeName(),
        *this,
        r,
        IOobject::NO_WRITE,
        false
    );
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldIOobject
(
    const Cloud<ParticleType>& c,
    const IOField<DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldFieldIOobject
(
    const Cloud<ParticleType>& c,
    const CompactIOField<Field<DataType>, DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeFields() const
{
    ParticleType::writeFields(*this);
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool
) const
{
    writeCloudUniformProperties();

    writeFields();
    return cloud::writeObject(fmt, ver, cmp, this->size());
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Cloud<ParticleType>& pc)
{
    pc.writeData(os);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Cloud<ParticleType>&)");

    return os;
}


// ************************************************************************* //
