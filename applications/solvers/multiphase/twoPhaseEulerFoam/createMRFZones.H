    IOMRFZoneList MRF(mesh);
    MRF.correctBoundaryVelocity(U1);
    MRF.correctBoundaryVelocity(U2);
    MRF.correctBoundaryVelocity(U);
