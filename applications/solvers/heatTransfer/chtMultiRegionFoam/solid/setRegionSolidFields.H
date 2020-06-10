solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

tmp<volScalarField> tCv = thermo.Cv();
const volScalarField& Cv = tCv();

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCv = thermo.Kappa()/Cv;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.R().transformVector(tkappaByCv());
    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& e = thermo.he();

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];

solidNoLoopControl& pimple = pimples.solid(i);
