surfaceScalarField& phi = fluid.phi();

const IOMRFZoneList& MRF = fluid.MRF();
Foam::fvModels& fvModels = fluid.fvModels(mesh);
Foam::fvConstraints& fvConstraints = fluid.fvConstraints(mesh);
