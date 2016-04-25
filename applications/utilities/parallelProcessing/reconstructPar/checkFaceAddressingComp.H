{
    // Foam version 2.1 changes the addressing of faces in faceProcAddressing
    // The following code checks and modifies the addressing for cases where
    // the decomposition has been done with the foam2.0 and earlier tools, but
    // the reconstruction is attempted with version 2.1 or later

    label minFaceIndex = labelMax;

    PtrList<labelIOList>& faceProcAddressing = procMeshes.faceProcAddressing();

    forAll(faceProcAddressing, proci)
    {
        const labelList& curFaceAddr = faceProcAddressing[proci];

        forAll(curFaceAddr, facei)
        {
            if (mag(curFaceAddr[facei]) < minFaceIndex)
            {
                minFaceIndex = mag(curFaceAddr[facei]);
            }
        }
    }

    if (minFaceIndex < 1)
    {
        WarningInFunction
            << "parallel decomposition addressing." << endl
            << "It looks like you are trying to reconstruct the case "
            << "decomposed with an earlier version of FOAM, which could\n"
            << "potentially cause compatibility problems.  The code will "
            << "attempt to update the addressing automatically; in case of\n"
            << "failure, please repeat the decomposition of the case using "
            << "the current version fo decomposePar"
            << endl;

        forAll(faceProcAddressing, proci)
        {
            labelList& curFaceAddr = faceProcAddressing[proci];

            forAll(curFaceAddr, facei)
            {
                curFaceAddr[facei] += sign(curFaceAddr[facei]);
            }

            faceProcAddressing[proci].write();
        }
    }
}
