    IOMRFZoneList mrfZones(mesh);

    forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
    {
        mrfZones.correctBoundaryVelocity(iter().U());
    }

    mrfZones.correctBoundaryVelocity(U);
