// write time values to case file

scalar timeCorrection = 0;
if (timeDirs[0].value() < 0)
{
    timeCorrection = - timeDirs[0].value();
    Info<< "Correcting time values. Adding " << timeCorrection << endl;
}

// the case file is always ASCII
Info << "write case: " << caseFileName.c_str() << endl;

OFstream caseFile(ensightDir/caseFileName, IOstream::ASCII);
caseFile.setf(ios_base::left);

caseFile.setf(ios_base::scientific, ios_base::floatfield);
caseFile.precision(5);

caseFile
    << "FORMAT" << nl
    << setw(16) << "type:" << "ensight gold" << nl << nl;

if (hasMovingMesh)
{
    caseFile
        << "GEOMETRY" << nl
        << setw(16) << "model: 1" << (dataMask/geometryName).c_str() << nl;
}
else
{
    caseFile
        << "GEOMETRY" << nl
        << setw(16) << "model:" << geometryName << nl;
}


// add information for clouds
// multiple clouds currently require the same time index
forAllConstIter(HashTable<HashTable<word>>, cloudFields, cloudIter)
{
    const word& cloudName = cloudIter.key();

    caseFile
        << setw(16) << "measured: 2"
        << fileName(dataMask/cloud::prefix/cloudName/"positions").c_str()
            << nl;
}
caseFile
    << nl << "VARIABLE" << nl;

forAllConstIter(HashTable<word>, volumeFields, fieldIter)
{
    const word& fieldName = fieldIter.key();
    const word& fieldType = fieldIter();
    string ensightType;

    if (fieldType == volScalarField::typeName)
    {
        ensightType = "scalar";
    }
    else if (fieldType == volVectorField::typeName)
    {
        ensightType = "vector";
    }
    else if (fieldType == volSphericalTensorField::typeName)
    {
        ensightType = "tensor symm";
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        ensightType = "tensor symm";
    }
    else if (fieldType == volTensorField::typeName)
    {
        ensightType = "tensor asym";
    }

    if (ensightType.size())
    {
        ensightCaseEntry
        (
            caseFile,
            ensightType,
            fieldName,
            dataMask
        );
    }
}

// TODO: allow similar/different time-steps for each cloud


label cloudNo = 0;
forAllConstIter(HashTable<HashTable<word>>, cloudFields, cloudIter)
{
    const word& cloudName = cloudIter.key();

    forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
    {
        const word& fieldName = fieldIter.key();
        const word& fieldType = fieldIter();
        string ensightType;

        if (fieldType == scalarIOField::typeName)
        {
            ensightType = "scalar";
        }
        else if (fieldType == vectorIOField::typeName)
        {
            ensightType = "vector";
        }
        else if (fieldType == tensorIOField::typeName)
        {
            ensightType = "tensor";
        }

        if (ensightType.size())
        {
            ensightCaseEntry
            (
                caseFile,
                ensightType,
                fieldName,
                dataMask,
                cloud::prefix/cloudName,
                cloudNo,
                2
            );
        }
    }
    cloudNo++;
}


// add time values
caseFile << nl << "TIME" << nl;

// time set 1 - geometry and volume fields
if (fieldTimesUsed.size())
{
    caseFile
        << "time set:        " << 1 << nl
        << "number of steps: " << fieldTimesUsed.size() << nl
        << "filename numbers:" << nl;

    label count = 0;
    forAll(fieldTimesUsed, i)
    {
        caseFile
            << " " << setw(12) << fieldTimesUsed[i];

        if (++count % 6 == 0)
        {
            caseFile << nl;
        }
    }

    caseFile
        << nl << "time values:" << nl;

    count = 0;
    forAll(fieldTimesUsed, i)
    {
        caseFile
            << " " << setw(12)
                << timeIndices[fieldTimesUsed[i]] + timeCorrection;

        if (++count % 6 == 0)
        {
            caseFile << nl;
        }
    }
    caseFile << nl << nl;
}

// TODO: allow similar/different time-steps for each cloud
cloudNo = 0;
forAllConstIter(HashTable<DynamicList<label>>, cloudTimesUsed, cloudIter)
{
    // const word& cloudName = cloudIter.key();
    const DynamicList<label>& timesUsed = cloudIter();

    if (timesUsed.size() && cloudNo == 0)
    {
        caseFile
            << "time set:        " << 2 << nl
            << "number of steps: " << timesUsed.size() << nl
            << "filename numbers:" << nl;

        label count = 0;
        forAll(timesUsed, i)
        {
            caseFile
                << " " << setw(12) << timesUsed[i];

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }

        caseFile
            << nl << "time values:" << nl;

        count = 0;
        forAll(timesUsed, i)
        {
            caseFile
                << " " << setw(12)
                    << timeIndices[timesUsed[i]] + timeCorrection;

            if (++count % 6 == 0)
            {
                caseFile << nl;
            }
        }
        caseFile << nl << nl;

        cloudNo++;
    }
}

caseFile << "# end" << nl;
