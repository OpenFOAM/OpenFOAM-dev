// Read surface fields.

PtrList<surfaceScalarField> ssFlds;
if (fields) ReadFields(mesh, objects, ssFlds);

PtrList<surfaceVectorField> svFlds;
if (fields) ReadFields(mesh, objects, svFlds);

PtrList<surfaceSphericalTensorField> sstFlds;
if (fields) ReadFields(mesh, objects, sstFlds);

PtrList<surfaceSymmTensorField> ssymtFlds;
if (fields) ReadFields(mesh, objects, ssymtFlds);

PtrList<surfaceTensorField> stFlds;
if (fields) ReadFields(mesh, objects, stFlds);
