{
    // Evaluate near-wall behaviour

    const scalar nu = turbulence->nu()().boundaryField()[patchId][faceId];
    const scalar nut = turbulence->nut()().boundaryField()[patchId][faceId];
    const scalar nuEff = turbulence->nuEff()().boundaryField()[patchId][faceId];
    const symmTensor R = turbulence->sigma()().boundaryField()[patchId][faceId];
    const scalar epsilon = turbulence->epsilon()()[cellId];
    const scalar omega = turbulence->omega()()[cellId];
    const scalar k = turbulence->k()()[cellId];
    const scalar magUp = mag(U[cellId] - U.boundaryField()[patchId][faceId]);

    const scalar tauw = flowDirection & R & wallNormal;
    const scalar uTau = ::sqrt(mag(tauw));

    const scalar yPlus = uTau*y[cellId]/(nu + rootVSmall);
    const scalar uPlus = magUp/(uTau + rootVSmall);
    const scalar nutPlus = nut/nu;
    const scalar nuEffPlus = nuEff/nu;
    const scalar kPlus = k/(sqr(uTau) + rootVSmall);
    const scalar epsilonPlus = epsilon*nu/(pow4(uTau) + rootVSmall);
    const scalar omegaPlus = omega*nu/(sqr(uTau) + rootVSmall);

    const scalar Rey = magUp*y[cellId]/nu;

    Info<< "Rey = " << Rey
        << ", uTau = " << uTau
        << ", nut+ = " << nutPlus
        << ", nuEff+ = " << nuEffPlus
        << ", y+ = " << yPlus << ", u+ = " << uPlus
        << ", k+ = " << kPlus << ", epsilon+ = " << epsilonPlus
        << ", omega+ = " << omegaPlus
        << endl;
}
