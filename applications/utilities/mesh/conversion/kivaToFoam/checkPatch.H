    if (bci != 0)
    {
        label bcIndex = bcIDs[bci];
        label regionIndex = 0;

        // Decompose moving faces into piston and valves using the idface array
        if (bci == 1)
        {
            if (kivaVersion == kiva3v)
            {
                regionIndex = max
                (
                    max(idface[quadFace[0]], idface[quadFace[1]]),
                    max(idface[quadFace[2]], idface[quadFace[3]])
                );

                if (regionIndex > 0)
                {
                    bcIndex = VALVE;
                    regionIndex--;
                }
            }
        }

        // Decompose fixed wall faces into cylinder-head and liner using
        // the fv array and split of the y=0 faces as a symmetry plane
        if (bci == 2)
        {
            if
            (
                fv[quadFace[0]] == 6
             || fv[quadFace[1]] == 6
             || fv[quadFace[2]] == 6
             || fv[quadFace[3]] == 6
            )
            {
                bcIndex = CYLINDERHEAD;
            }

            if
            (
                   mag(points[quadFace[0]].y()) < small
                && mag(points[quadFace[1]].y()) < small
                && mag(points[quadFace[2]].y()) < small
                && mag(points[quadFace[3]].y()) < small
            )
            {
                bcIndex = SYMMETRYPLANE;
            }
        }

        // Make the back part of the periodic boundary the second region
        // of the front part
        if (bci == 6)
        {
            bcIndex = WEDGE;
            regionIndex = 1;
        }

        if (regionIndex >= pFaces[bcIndex].size())
        {
            pFaces[bcIndex].setSize(regionIndex + 1);
        }

        quadFace[0] = pointMap[quadFace[0]];
        quadFace[1] = pointMap[quadFace[1]];
        quadFace[2] = pointMap[quadFace[2]];
        quadFace[3] = pointMap[quadFace[3]];

        label ti = 0;
        forAll(quadFace, i)
        {
            if (quadFace[i] != quadFace[(i+1)%4])
            {
                triFace[(ti++)%3] = quadFace[i];
            }
        }

        if (ti == 4)
        {
            pFaces[bcIndex][regionIndex].append(quadFace);
        }
        else if (ti == 3)
        {
            pFaces[bcIndex][regionIndex].append(triFace);
        }
        // else the face has collapsed to an edge or point
    }
    else
    {
        FatalErrorInFunction
            << "bc not defined for active cell = " << i
            << exit(FatalError);
    }
