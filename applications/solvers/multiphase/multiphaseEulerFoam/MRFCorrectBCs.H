    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        MRF.correctBoundaryVelocity(iter().U());
    }

    MRF.correctBoundaryVelocity(U);
