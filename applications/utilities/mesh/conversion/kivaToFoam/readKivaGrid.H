ifstream kivaFile(kivaFileName.c_str());

if (!kivaFile.good())
{
    FatalErrorInFunction
        << "Cannot open file " << kivaFileName
        << exit(FatalError);
}

Info << "Reading kiva grid from file " << kivaFileName << endl;


char title[120];
kivaFile.getline(title, 120, '\n');

label nPoints, nCells, nRegs;
kivaFile >> nCells >> nPoints >> nRegs;

pointField points(nPoints);
label i4;
labelList idface(nPoints), fv(nPoints);

for (label i=0; i<nPoints; i++)
{
    scalar ffv;
    kivaFile
        >> i4
        >> points[i].x() >> points[i].y() >> points[i].z()
        >> ffv;

    if (kivaVersion == kiva3v)
    {
        kivaFile >> idface[i];
    }

    // Convert from KIVA cgs units to SI
    points[i] *= 0.01;

    fv[i] = label(ffv);
}

labelList i1tab(nPoints), i3tab(nPoints), i8tab(nPoints), idreg(nPoints),
      f(nPoints), bcl(nPoints), bcf(nPoints), bcb(nPoints);

label nBfaces = 0;

for (label i=0; i<nPoints; i++)
{
    label i1, i3, i8;
    scalar ff, fbcl, fbcf, fbcb;

    kivaFile
        >> i1 >> i3 >> i4 >> i8
        >> ff >> fbcl >> fbcf >> fbcb
        >> idreg[i];

    // Correct indices to start from 0
    i1tab[i] = i1 - 1;
    i3tab[i] = i3 - 1;
    i8tab[i] = i8 - 1;

    // Convert scalar indices into hexLabels
    f[i] = label(ff);
    bcl[i] = label(fbcl);
    bcf[i] = label(fbcf);
    bcb[i] = label(fbcb);

    if (bcl[i] > 0 && bcl[i] != 4) nBfaces++;
    if (bcf[i] > 0 && bcf[i] != 4) nBfaces++;
    if (bcb[i] > 0 && bcb[i] != 4) nBfaces++;
}


label mTable;
kivaFile >> mTable;

if (mTable == 0)
{
    FatalErrorInFunction
        << "mTable == 0. This is not supported."
           " kivaToFoam needs complete neighbour information"
        << exit(FatalError);
}


labelList imtab(nPoints), jmtab(nPoints), kmtab(nPoints);

for (label i=0; i<nPoints; i++)
{
    label im, jm, km;
    kivaFile >> i4 >> im >> jm >> km;

    // Correct indices to start from 0
    imtab[i] = im - 1;
    jmtab[i] = jm - 1;
    kmtab[i] = km - 1;
}

Info << "Finished reading KIVA file" << endl;

cellShapeList cellShapes(nPoints);

labelList cellZoning(nPoints, -1);

const cellModel& hex = *(cellModeller::lookup("hex"));
labelList hexLabels(8);
label activeCells = 0;

// Create and set the collocated point collapse map
labelList pointMap(nPoints);

forAll(pointMap, i)
{
    pointMap[i] = i;
}

// Initialise all cells to hex and search and set map for collocated points
for (label i=0; i<nPoints; i++)
{
    if (f[i] > 0.0)
    {
        hexLabels[0] = i;
        hexLabels[1] = i1tab[i];
        hexLabels[2] = i3tab[i1tab[i]];
        hexLabels[3] = i3tab[i];
        hexLabels[4] = i8tab[i];
        hexLabels[5] = i1tab[i8tab[i]];
        hexLabels[6] = i3tab[i1tab[i8tab[i]]];
        hexLabels[7] = i3tab[i8tab[i]];

        cellShapes[activeCells] = cellShape(hex, hexLabels);

        edgeList edges = cellShapes[activeCells].edges();

        forAll(edges, ei)
        {
            if (edges[ei].mag(points) < small)
            {
                label start = pointMap[edges[ei].start()];
                while (start != pointMap[start])
                {
                    start = pointMap[start];
                }

                label end = pointMap[edges[ei].end()];
                while (end != pointMap[end])
                {
                    end = pointMap[end];
                }

                label minLabel = min(start, end);

                pointMap[start] = pointMap[end] = minLabel;
            }
        }

        // Grab the cell zoning info
        cellZoning[activeCells] = idreg[i];

        activeCells++;
    }
}

// Contract list of cells to the active ones
cellShapes.setSize(activeCells);
cellZoning.setSize(activeCells);

// Map collocated points to refer to the same point and collapse cell shape
// to the corresponding hex-degenerate.
forAll(cellShapes, celli)
{
    cellShape& cs = cellShapes[celli];

    forAll(cs, i)
    {
        cs[i] = pointMap[cs[i]];
    }

    cs.collapse();
}

label bcIDs[11] = {-1, 0, 2, 4, -1, 5, -1, 6, 7, 8, 9};

const label nBCs = 12;

const word* kivaPatchTypes[nBCs] =
{
    &wallPolyPatch::typeName,
    &wallPolyPatch::typeName,
    &wallPolyPatch::typeName,
    &wallPolyPatch::typeName,
    &symmetryPolyPatch::typeName,
    &wedgePolyPatch::typeName,
    &polyPatch::typeName,
    &polyPatch::typeName,
    &polyPatch::typeName,
    &polyPatch::typeName,
    &symmetryPolyPatch::typeName,
    &mergedCyclicPolyPatch::typeName
};

enum patchTypeNames
{
    PISTON,
    VALVE,
    LINER,
    CYLINDERHEAD,
    AXIS,
    WEDGE,
    INFLOW,
    OUTFLOW,
    PRESIN,
    PRESOUT,
    SYMMETRYPLANE,
    CYCLIC
};

const char* kivaPatchNames[nBCs] =
{
    "piston",
    "valve",
    "liner",
    "cylinderHead",
    "axis",
    "wedge",
    "inflow",
    "outflow",
    "presin",
    "presout",
    "symmetryPlane",
    "cyclic"
};


List<SLList<face>> pFaces[nBCs];

face quadFace(4);
face triFace(3);

for (label i=0; i<nPoints; i++)
{
    if (f[i] > 0)
    {
        // left
        label bci = bcl[i];
        if (bci != 4)
        {
            quadFace[0] = i3tab[i];
            quadFace[1] = i;
            quadFace[2] = i8tab[i];
            quadFace[3] = i3tab[i8tab[i]];

            #include "checkPatch.H"
        }

        // right
        bci = bcl[i1tab[i]];
        if (bci != 4)
        {
            quadFace[0] = i1tab[i];
            quadFace[1] = i3tab[i1tab[i]];
            quadFace[2] = i8tab[i3tab[i1tab[i]]];
            quadFace[3] = i8tab[i1tab[i]];

            #include "checkPatch.H"
        }

        // front
        bci = bcf[i];
        if (bci != 4)
        {
            quadFace[0] = i;
            quadFace[1] = i1tab[i];
            quadFace[2] = i8tab[i1tab[i]];
            quadFace[3] = i8tab[i];

            #include "checkPatch.H"
        }

        // derriere (back)
        bci = bcf[i3tab[i]];
        if (bci != 4)
        {
            quadFace[0] = i3tab[i];
            quadFace[1] = i8tab[i3tab[i]];
            quadFace[2] = i8tab[i1tab[i3tab[i]]];
            quadFace[3] = i1tab[i3tab[i]];

            #include "checkPatch.H"
        }

        // bottom
        bci = bcb[i];
        if (bci != 4)
        {
            quadFace[0] = i;
            quadFace[1] = i3tab[i];
            quadFace[2] = i1tab[i3tab[i]];
            quadFace[3] = i1tab[i];

            #include "checkPatch.H"
        }

        // top
        bci = bcb[i8tab[i]];
        if (bci != 4)
        {
            quadFace[0] = i8tab[i];
            quadFace[1] = i1tab[i8tab[i]];
            quadFace[2] = i3tab[i1tab[i8tab[i]]];
            quadFace[3] = i3tab[i8tab[i]];

            #include "checkPatch.H"
        }
    }
}

// Transfer liner faces that are above the minimum cylinder-head z height
// into the cylinder-head region
if
(
    pFaces[LINER].size()
 && pFaces[LINER][0].size()
 && pFaces[CYLINDERHEAD].size()
 && pFaces[CYLINDERHEAD][0].size()
)
{
    scalar minz = great;

    forAllConstIter(SLList<face>, pFaces[CYLINDERHEAD][0], iter)
    {
        const face& pf = iter();

        forAll(pf, pfi)
        {
            minz = min(minz, points[pf[pfi]].z());
        }
    }

    minz += small;

    SLList<face> newLinerFaces;

    forAllConstIter(SLList<face>, pFaces[LINER][0], iter)
    {
        const face& pf = iter();

        scalar minfz = great;
        forAll(pf, pfi)
        {
            minfz = min(minfz, points[pf[pfi]].z());
        }

        if (minfz > minz)
        {
            pFaces[CYLINDERHEAD][0].append(pf);
        }
        else
        {
            newLinerFaces.append(pf);
        }
    }

    if (pFaces[LINER][0].size() != newLinerFaces.size())
    {
        Info<< "Transferred " << pFaces[LINER][0].size() - newLinerFaces.size()
            << " faces from liner region to cylinder head" << endl;
        pFaces[LINER][0] = newLinerFaces;
    }

    SLList<face> newCylinderHeadFaces;

    forAllConstIter(SLList<face>, pFaces[CYLINDERHEAD][0], iter)
    {
        const face& pf = iter();

        scalar minfz = great;
        forAll(pf, pfi)
        {
            minfz = min(minfz, points[pf[pfi]].z());
        }

        if (minfz < zHeadMin)
        {
            pFaces[LINER][0].append(pf);
        }
        else
        {
            newCylinderHeadFaces.append(pf);
        }
    }

    if (pFaces[CYLINDERHEAD][0].size() != newCylinderHeadFaces.size())
    {
        Info<< "Transferred faces from cylinder-head region to linder" << endl;
        pFaces[CYLINDERHEAD][0] = newCylinderHeadFaces;
    }
}


// Count the number of non-zero sized patches
label nPatches = 0;
for (int bci=0; bci<nBCs; bci++)
{
    forAll(pFaces[bci], rgi)
    {
        if (pFaces[bci][rgi].size())
        {
            nPatches++;
        }
    }
}


// Sort-out the 2D/3D wedge geometry
if (pFaces[WEDGE].size() && pFaces[WEDGE][0].size())
{
    if (pFaces[WEDGE][0].size() == cellShapes.size())
    {
        // Distribute the points to be +/- 2.5deg from the x-z plane

        scalar tanTheta = Foam::tan(degToRad(2.5));

        SLList<face>::iterator iterf = pFaces[WEDGE][0].begin();
        SLList<face>::iterator iterb = pFaces[WEDGE][1].begin();
        for
        (
            ;
            iterf != pFaces[WEDGE][0].end() && iterb != pFaces[WEDGE][1].end();
            ++iterf, ++iterb
        )
        {
            for (direction d=0; d<4; d++)
            {
                points[iterf()[d]].y() = -tanTheta*points[iterf()[d]].x();
                points[iterb()[d]].y() =  tanTheta*points[iterb()[d]].x();
            }
        }
    }
    else
    {
        pFaces[CYCLIC].setSize(1);
        pFaces[CYCLIC][0] = pFaces[WEDGE][0];
        forAllIter(SLList<face>, pFaces[WEDGE][1], iterb)
        {
            pFaces[CYCLIC][0].append(iterb());
        }

        pFaces[WEDGE].clear();
        nPatches--;
    }
}


// Build the patches

faceListList boundary(nPatches);
wordList patchNames(nPatches);
wordList patchTypes(nPatches);
word defaultFacesName = "defaultFaces";
word defaultFacesType = emptyPolyPatch::typeName;

label nAddedPatches = 0;

for (int bci=0; bci<nBCs; bci++)
{
    forAll(pFaces[bci], rgi)
    {
        if (pFaces[bci][rgi].size())
        {
            patchTypes[nAddedPatches] = *kivaPatchTypes[bci];

            patchNames[nAddedPatches] = kivaPatchNames[bci];

            if (pFaces[bci].size() > 1)
            {
                patchNames[nAddedPatches] += name(rgi+1);
            }

            boundary[nAddedPatches] = pFaces[bci][rgi];
            nAddedPatches++;
        }
    }
}


// Remove unused points and update cell and face labels accordingly

labelList pointLabels(nPoints, -1);

// Scan cells for used points
forAll(cellShapes, celli)
{
    forAll(cellShapes[celli], i)
    {
        pointLabels[cellShapes[celli][i]] = 1;
    }
}

// Create addressing for used points and pack points array
label newPointi = 0;
forAll(pointLabels, pointi)
{
    if (pointLabels[pointi] != -1)
    {
        pointLabels[pointi] = newPointi;
        points[newPointi++] = points[pointi];
    }
}
points.setSize(newPointi);

// Reset cell point labels
forAll(cellShapes, celli)
{
    cellShape& cs = cellShapes[celli];

    forAll(cs, i)
    {
        cs[i] = pointLabels[cs[i]];
    }
}

// Reset boundary-face point labels
forAll(boundary, patchi)
{
    forAll(boundary[patchi], facei)
    {
        face& f = boundary[patchi][facei];

        forAll(f, i)
        {
            f[i] = pointLabels[f[i]];
        }
    }
}

PtrList<dictionary> patchDicts;
preservePatchTypes
(
    runTime,
    runTime.constant(),
    polyMesh::meshSubDir,
    patchNames,
    patchDicts,
    defaultFacesName,
    defaultFacesType
);
// Add information to dictionary
forAll(patchNames, patchi)
{
    if (!patchDicts.set(patchi))
    {
        patchDicts.set(patchi, new dictionary());
    }
    // Add but not overwrite
    patchDicts[patchi].add("type", patchTypes[patchi], false);
}

// Build the mesh and write it out
polyMesh pShapeMesh
(
    IOobject
    (
        polyMesh::defaultRegion,
        runTime.constant(),
        runTime
    ),
    move(points),
    cellShapes,
    boundary,
    patchNames,
    patchDicts,
    defaultFacesName,
    defaultFacesType
);

// Un-merge any merged cyclics
polyMeshUnMergeCyclics(pShapeMesh);

Info << "Writing polyMesh" << endl;
pShapeMesh.write();

fileName czPath
(
    runTime.path()/runTime.constant()/polyMesh::defaultRegion/"cellZoning"
);
Info << "Writing cell zoning info to file: " << czPath << endl;
OFstream cz(czPath);

cz << cellZoning << endl;
