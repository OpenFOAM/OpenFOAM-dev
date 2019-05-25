Info<< "Constructing single cell mesh" << nl << endl;

labelList owner(6, label(0));
labelList neighbour(0);

pointField points(8);
points[0] = vector(0, 0, 0);
points[1] = vector(1, 0, 0);
points[2] = vector(1, 1, 0);
points[3] = vector(0, 1, 0);
points[4] = vector(0, 0, 1);
points[5] = vector(1, 0, 1);
points[6] = vector(1, 1, 1);
points[7] = vector(0, 1, 1);

const cellModel& hexa = *(cellModeller::lookup("hex"));
faceList faces = hexa.modelFaces();

fvMesh mesh
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::READ_IF_PRESENT
    ),
    move(points),
    move(faces),
    move(owner),
    move(neighbour)
);

List<polyPatch*> patches(1);

patches[0] = new emptyPolyPatch
(
    "boundary",
    6,
    0,
    0,
    mesh.boundaryMesh(),
    emptyPolyPatch::typeName
);

mesh.addFvPatches(patches);
