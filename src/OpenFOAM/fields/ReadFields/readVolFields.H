// Read vol fields.

PtrList<volScalarField> vsFlds;
if (fields) ReadFields(mesh, objects, vsFlds);

PtrList<volVectorField> vvFlds;
if (fields) ReadFields(mesh, objects, vvFlds);

PtrList<volSphericalTensorField> vstFlds;
if (fields) ReadFields(mesh, objects, vstFlds);

PtrList<volSymmTensorField> vsymtFlds;
if (fields) ReadFields(mesh, objects, vsymtFlds);

PtrList<volTensorField> vtFlds;
if (fields) ReadFields(mesh, objects, vtFlds);
