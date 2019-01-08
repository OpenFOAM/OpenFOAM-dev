// Read point fields.

PtrList<pointScalarField> psFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, psFlds);

PtrList<pointVectorField> pvFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, pvFlds);

PtrList<pointSphericalTensorField> pstFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, pstFlds);

PtrList<pointSymmTensorField> psymtFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, psymtFlds);

PtrList<pointTensorField> ptFlds;
if (fields) ReadFields(pointMesh::New(mesh), objects, ptFlds);
