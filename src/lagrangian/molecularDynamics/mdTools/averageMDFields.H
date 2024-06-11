if (runTime.writeTime())
{
    /*-----------------------------------------------------------------------*\
        Number density
    \*-----------------------------------------------------------------------*/

    scalarField totalRhoN_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesRhoN, rN)
    {
        allSpeciesRhoN[rN].primitiveFieldRef() =
            allSpeciesN_RU[rN]
            /mesh.cellVolumes()
            /nAveragingSteps;

        totalRhoN_sum += allSpeciesRhoN[rN].primitiveField();
    }

    totalRhoN.primitiveFieldRef() = totalRhoN_sum;


    /*-----------------------------------------------------------------------*\
        Mass density
    \*-----------------------------------------------------------------------*/

    scalarField totalRhoM_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesRhoM, rM)
    {
        allSpeciesRhoM[rM].primitiveFieldRef() =
            allSpeciesM_RU[rM]
            /mesh.cellVolumes()
            /nAveragingSteps;

        totalRhoM_sum += allSpeciesRhoM[rM].primitiveField();
    }

    totalRhoM.primitiveFieldRef() = totalRhoM_sum;

    /*-----------------------------------------------------------------------*\
        Bulk velocity
    \*-----------------------------------------------------------------------*/

    vectorField totalMomentum_sum(mesh.nCells(), Zero);

    scalarField totalMass_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesVelocity, v)
    {
        // A check for 1/0 molecules is required.

        vectorField& singleSpeciesVelocity
        (
            allSpeciesVelocity[v].primitiveField()
        );

        forAll(singleSpeciesVelocity, sSV)
        {
            if (allSpeciesN_RU[v][sSV])
            {
                singleSpeciesVelocity[sSV] =
                    allSpeciesVelocitySum_RU[v][sSV]
                    /allSpeciesN_RU[v][sSV];

                totalMomentum_sum[sSV] +=
                    allSpeciesM_RU[v][sSV]
                    /allSpeciesN_RU[v][sSV]
                    *allSpeciesVelocitySum_RU[v][sSV];

                totalMass_sum[sSV] += allSpeciesM_RU[v][sSV];
            }
            else
            {
                singleSpeciesVelocity[sSV] = Zero;
            }
        }
    }

    volVectorField::Internal& itotalVelocity = totalVelocity;

    forAll(itotalVelocity, tV)
    {
        if (totalMass_sum[tV] > vSmall)
        {
            itotalVelocity[tV] = totalMomentum_sum[tV]/totalMass_sum[tV];
        }
        else
        {
            itotalVelocity[tV] = Zero;
        }
    }

    /*-----------------------------------------------------------------------*\
        Kinetic temperature
    \*-----------------------------------------------------------------------*/

    scalarField totalTemperatureVTerms_sum(mesh.nCells(), 0.0);

    scalarField totalN_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesTemperature, t)
    {
        // A check for 1/0 molecules is required.

        scalarField& singleSpeciesTemp
        (
            allSpeciesTemperature[t].primitiveField()
        );

        forAll(singleSpeciesTemp, sST)
        {
            if (allSpeciesN_RU[t][sST])
            {
                singleSpeciesTemp[sST] =
                    allSpeciesM_RU[t][sST]
                    /allSpeciesN_RU[t][sST]
                    /(3.0 * moleculeCloud::kb * allSpeciesN_RU[t][sST])
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[t][sST]
                        -
                        (
                            allSpeciesVelocitySum_RU[t][sST]
                            &
                            allSpeciesVelocitySum_RU[t][sST]
                        )
                        /allSpeciesN_RU[t][sST]
                    );

                totalTemperatureVTerms_sum[sST] +=
                    allSpeciesM_RU[t][sST]
                   /allSpeciesN_RU[t][sST]
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[t][sST]
                      -
                        (
                            allSpeciesVelocitySum_RU[t][sST]
                            &
                            allSpeciesVelocitySum_RU[t][sST]
                        )
                        /allSpeciesN_RU[t][sST]
                    );

                totalN_sum[sST] += allSpeciesN_RU[t][sST];
            }
            else
            {
                singleSpeciesTemp[sST] = 0.0;
            }
        }
    }

    volScalarField::Internal& itotalTemperature =
        totalTemperature;

    forAll(itotalTemperature, tT)
    {
        if (totalN_sum[tT] > 0)
        {
            itotalTemperature[tT] =
                totalTemperatureVTerms_sum[tT]
                /(3.0 * moleculeCloud::kb * totalN_sum[tT]);
        }
        else
        {
            itotalTemperature[tT] = 0.0;
        }
    }

    /*-----------------------------------------------------------------------*\
        Mean kinetic energy
    \*-----------------------------------------------------------------------*/

    scalarField totalKE_sum(mesh.nCells(), 0.0);

    forAll(allSpeciesMeanKE, mKE)
    {
        // A check for 1/0 molecules is required.

        scalarField& singleSpeciesMeanKE
        (
            allSpeciesMeanKE[mKE].primitiveField()
        );

        forAll(singleSpeciesMeanKE, sSMKE)
        {
            if (allSpeciesN_RU[mKE][sSMKE])
            {
                singleSpeciesMeanKE[sSMKE] =
                    allSpeciesM_RU[mKE][sSMKE]
                   /allSpeciesN_RU[mKE][sSMKE]
                   /(2.0*allSpeciesN_RU[mKE][sSMKE])
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[mKE][sSMKE]
                    );

                totalKE_sum[sSMKE] +=
                    allSpeciesM_RU[mKE][sSMKE]
                    /allSpeciesN_RU[mKE][sSMKE]
                    /2.0
                   *(
                        allSpeciesVelocityMagSquaredSum_RU[mKE][sSMKE]
                    );
            }
            else
            {
                singleSpeciesMeanKE[sSMKE] = 0.0;
            }
        }
    }

    volScalarField::Internal& itotalMeanKE = totalMeanKE;

    forAll(itotalMeanKE, tMKE)
    {
        if (totalN_sum[tMKE] > 0)
        {
            itotalMeanKE[tMKE] =
                totalKE_sum[tMKE]
                /totalN_sum[tMKE];
        }
        else
        {
            itotalMeanKE[tMKE] = 0.0;
        }
    }
}
