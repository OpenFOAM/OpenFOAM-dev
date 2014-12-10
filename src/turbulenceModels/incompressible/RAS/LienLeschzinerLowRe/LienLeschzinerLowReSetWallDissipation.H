// Set the dissipation in the near-wall cell to the value prescribed by the
// Lien-Leschziner low-Re model with Wolfstein length-scale prescription

{
    labelList cellBoundaryFaceCount(epsilon_.size(), 0);

    //- use constant Cmu for epsilon in the near-wall cell
    scalar Cmu75 = pow(Cmu_.value(), 0.75);

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall epsilon field to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] = 0.0;
            }
        }
    }

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                // For corner cells (with two boundary or more faces),
                // epsilon in the near-wall cell are calculated as an average

                cellBoundaryFaceCount[faceCelli]++;

                epsilon_[faceCelli] +=
                     Cmu75*pow(k_[faceCelli], 1.5)
                    /(
                         kappa_.value()*y_[faceCelli]
                         *(1.0 - exp(-Aepsilon_.value()*yStar_[faceCelli]))
                     )
                    *exp(-Amu_.value()*sqr(yStar_[faceCelli]));

            }
        }
    }

    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}
