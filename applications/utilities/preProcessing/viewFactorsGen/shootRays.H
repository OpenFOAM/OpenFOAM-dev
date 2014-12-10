// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.

// Maximum lenght for dynamicList
const label maxDynListLength = 100000;

for (label procI = 0; procI < Pstream::nProcs(); procI++)
{
    // Shoot rays from me to procI. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    DynamicField<point> start(coarseMesh.nFaces());
    DynamicField<point> end(start.size());
    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());

    DynamicList<label> startAgg(start.size());
    DynamicList<label> endAgg(start.size());


    const pointField& myFc = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];
    const labelField& myAgg = remoteCoarseAgg[Pstream::myProcNo()];

    const pointField& remoteArea = remoteCoarseSf[procI];
    const pointField& remoteFc = remoteCoarseCf[procI];
    const labelField& remoteAgg = remoteCoarseAgg[procI];

    label i = 0;
    label j = 0;
    do
    {
        for (; i < myFc.size(); i++)
        {
            const point& fc = myFc[i];
            const vector& fA = myArea[i];
            const label& fAgg = myAgg[i];

            for (; j < remoteFc.size(); j++)//
            {
                if (procI != Pstream::myProcNo() || i != j)
                {
                    const point& remFc = remoteFc[j];
                    const vector& remA = remoteArea[j];
                    const label& remAgg = remoteAgg[j];

                    const vector& d = remFc - fc;

                    if (((d & fA) < 0.) && ((d & remA) > 0))
                    {
                        start.append(fc + 0.001*d);
                        startIndex.append(i);
                        startAgg.append(globalNumbering.toGlobal(procI, fAgg));
                        end.append(fc + 0.999*d);
                        label globalI = globalNumbering.toGlobal(procI, j);
                        endIndex.append(globalI);
                        endAgg.append(globalNumbering.toGlobal(procI, remAgg));
                        if (startIndex.size() > maxDynListLength)
                        {
                            FatalErrorIn
                            (
                                "shootRays"
                            )   << "Dynamic list need from capacity."
                                << "Actual size maxDynListLength : "
                                <<  maxDynListLength
                                << abort(FatalError);
                        }
                    }
                }
            }

            if (j == remoteFc.size())
            {
                j = 0;
            }
        }

    }while (returnReduce(i < myFc.size(), orOp<bool>()));

    List<pointIndexHit> hitInfo(startIndex.size());
    surfacesMesh.findLine(start, end, hitInfo);

    // Return hit global agglo index
    labelList aggHitIndex;
    surfacesMesh.getField(hitInfo, aggHitIndex);

    DynamicList<label> dRayIs;

    // Collect the rays which has not abstacle in bettween in rayStartFace
    // and rayEndFace. If the ray hit itself get stored in dRayIs
    forAll (hitInfo, rayI)
    {
        if (!hitInfo[rayI].hit())
        {
            rayStartFace.append(startIndex[rayI]);
            rayEndFace.append(endIndex[rayI]);
        }
        else if (aggHitIndex[rayI] == startAgg[rayI])
        {
            dRayIs.append(rayI);
        }
    }

    start.clear();


    // Continue rays which hit themself. If they hit the target
    // agglomeration are added to rayStartFace and rayEndFace

    bool firstLoop = true;
    DynamicField<point> startHitItself;
    DynamicField<point> endHitItself;
    label iter = 0;

    do
    {
        labelField rayIs;
        rayIs.transfer(dRayIs);
        dRayIs.clear();
        forAll (rayIs, rayI)
        {
            const label rayID = rayIs[rayI];
            label hitIndex = -1;

            if (firstLoop)
            {
                hitIndex = rayIs[rayI];
            }
            else
            {
                hitIndex = rayI;
            }

            if (hitInfo[hitIndex].hit())
            {
                if (aggHitIndex[hitIndex] == startAgg[rayID])
                {
                    const vector& endP = end[rayID];
                    const vector& startP = hitInfo[hitIndex].hitPoint();
                    const vector& d = endP - startP;

                    startHitItself.append(startP + 0.01*d);
                    endHitItself.append(startP + 1.01*d);

                    dRayIs.append(rayID);
                }
                else if (aggHitIndex[hitIndex] == endAgg[rayID])
                {
                    rayStartFace.append(startIndex[rayID]);
                    rayEndFace.append(endIndex[rayID]);
                }

            }
        }

        hitInfo.clear();
        hitInfo.resize(dRayIs.size());

        surfacesMesh.findLine(startHitItself, endHitItself, hitInfo);

        surfacesMesh.getField(hitInfo, aggHitIndex);


        endHitItself.clear();
        startHitItself.clear();
        firstLoop = false;
        iter ++;

    }while (returnReduce(hitInfo.size(), orOp<bool>()) > 0 && iter < 10);

    startIndex.clear();
    end.clear();
    endIndex.clear();
    startAgg.clear();
    endAgg.clear();
}
