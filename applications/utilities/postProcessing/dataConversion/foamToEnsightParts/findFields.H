// check the final time directory for the following:

// 1. volume fields
HashTable<word> volumeFields;

// 2. the fields for each cloud:
HashTable< HashTable<word> > cloudFields;

if (timeDirs.size())
{
    IOobjectList objs(mesh, timeDirs.last().name());

    forAllConstIter(IOobjectList, objs, fieldIter)
    {
        const IOobject& obj = *fieldIter();
        const word& fieldName = obj.name();
        const word& fieldType = obj.headerClassName();

        if (fieldName.size() > 2 && fieldName(fieldName.size()-2, 2) == "_0")
        {
            // ignore _0 fields
        }
        else if (volFieldTypes.found(fieldType))
        {
            // simply ignore types that we don't handle
            volumeFields.insert(fieldName, fieldType);
        }
    }


    //
    // now check for lagrangian/<cloudName>
    //
    fileNameList cloudDirs = readDir
    (
        runTime.path()
      / timeDirs.last().name()
      / regionPrefix
      / cloud::prefix,
        fileName::DIRECTORY
    );

    forAll(cloudDirs, cloudI)
    {
        const word& cloudName = cloudDirs[cloudI];

        // Create a new hash table for each cloud
        cloudFields.insert(cloudName, HashTable<word>());

        // Identify the new cloud within the hash table
        HashTable<HashTable<word> >::iterator cloudIter =
            cloudFields.find(cloudName);

        IOobjectList objs
        (
            mesh,
            timeDirs.last().name(),
            cloud::prefix/cloudName
        );

        bool hasPositions = false;
        forAllConstIter(IOobjectList, objs, fieldIter)
        {
            const IOobject obj = *fieldIter();
            const word& fieldName = obj.name();
            const word& fieldType = obj.headerClassName();

            if (fieldName == "positions")
            {
                hasPositions = true;
            }
            else if (cloudFieldTypes.found(fieldType))
            {
                // simply ignore types that we don't handle
                cloudIter().insert(fieldName, fieldType);
            }
        }

        // drop this cloud if it has no positions or is otherwise empty
        if (!hasPositions || cloudIter().empty())
        {
            Info<< "removing cloud " << cloudName << endl;
            cloudFields.erase(cloudIter);
        }
    }

    //
    // verify that the variable is present for all times
    //
    for (label i=0; volumeFields.size() && i < timeDirs.size(); ++i)
    {
        IOobjectList objs(mesh, timeDirs[i].name());

        forAllIter(HashTable<word>, volumeFields, fieldIter)
        {
            const word& fieldName = fieldIter.key();

            if (!objs.found(fieldName))
            {
                volumeFields.erase(fieldIter);
            }
        }
    }
}

