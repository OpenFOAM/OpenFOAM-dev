    IOMRFZoneList mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U1);
    mrfZones.correctBoundaryVelocity(U2);
    mrfZones.correctBoundaryVelocity(U);
