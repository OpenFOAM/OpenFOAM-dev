{
    // Evaluate near-wall behaviour

    scalar nu = turbulence->nu()().boundaryField()[patchId][faceId];
    scalar nut = turbulence->nut()().boundaryField()[patchId][faceId];
    symmTensor R = turbulence->devReff()().boundaryField()[patchId][faceId];
    scalar epsilon = turbulence->epsilon()()[cellId];
//    scalar omega = turbulence->omega()()[cellId];
    scalar k = turbulence->k()()[cellId];
    scalar magUp = mag(U[cellId] - U.boundaryField()[patchId][faceId]);

    scalar tauw = flowDirection & R & wallNormal;

    scalar uTau = ::sqrt(mag(tauw));

    scalar yPlus = uTau*y[cellId]/(nu + ROOTVSMALL);

    scalar uPlus = magUp/(uTau + ROOTVSMALL);

    scalar nutPlus = nut/nu;

    scalar kPlus = k/(sqr(uTau) + ROOTVSMALL);

    scalar epsilonPlus = epsilon*nu/(pow4(uTau) + ROOTVSMALL);

//    scalar omegaPlus = omega*nu/(sqr(uTau) + ROOTVSMALL);

    scalar Rey = magUp*y[cellId]/nu;

    Info<< "Rey = " << Rey << ", uTau = " << uTau << ", nut+ = " << nutPlus
        << ", y+ = " << yPlus << ", u+ = " << uPlus
        << ", k+ = " << kPlus << ", epsilon+ = " << epsilonPlus
        << endl;
}
