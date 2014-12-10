/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "forces.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "fluidThermo.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forces, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::forces::createFileNames(const dictionary& dict) const
{
    DynamicList<word> names(1);

    const word forceType(dict.lookup("type"));

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        label nb = readLabel(binDict.lookup("nBin"));
        if (nb > 0)
        {
            names.append(forceType + "_bins");
        }
    }

    names.append(forceType);

    return names;
}


void Foam::forces::writeFileHeader(const label i)
{
    if (i == 0)
    {
        // force data

        writeHeader(file(i), "Forces");
        writeHeaderValue(file(i), "CofR", coordSys_.origin());
        writeCommented(file(i), "Time");

        file(i)
            << "forces(pressure viscous porous) "
            << "moment(pressure viscous porous)";

        if (localSystem_)
        {
            file(i)
                << tab
                << "localForces(pressure,viscous,porous) "
                << "localMoments(pressure,viscous,porous)";
        }
    }
    else if (i == 1)
    {
        // bin data

        writeHeader(file(i), "Force bins");
        writeHeaderValue(file(i), "bins", nBin_);
        writeHeaderValue(file(i), "start", binMin_);
        writeHeaderValue(file(i), "delta", binDx_);
        writeHeaderValue(file(i), "direction", binDir_);

        vectorField binPoints(nBin_);
        writeCommented(file(i), "x co-ords  :");
        forAll(binPoints, pointI)
        {
            binPoints[pointI] = (binMin_ + (pointI + 1)*binDx_)*binDir_;
            file(i) << tab << binPoints[pointI].x();
        }
        file(i) << nl;

        writeCommented(file(i), "y co-ords  :");
        forAll(binPoints, pointI)
        {
            file(i) << tab << binPoints[pointI].y();
        }
        file(i) << nl;

        writeCommented(file(i), "z co-ords  :");
        forAll(binPoints, pointI)
        {
            file(i) << tab << binPoints[pointI].z();
        }
        file(i) << nl;

        writeCommented(file(i), "Time");

        for (label j = 0; j < nBin_; j++)
        {
            const word jn('(' + Foam::name(j) + ')');
            const word f("forces" + jn + "[pressure,viscous,porous]");
            const word m("moments" + jn + "[pressure,viscous,porous]");

            file(i)<< tab << f << tab << m;
        }
        if (localSystem_)
        {
            for (label j = 0; j < nBin_; j++)
            {
                const word jn('(' + Foam::name(j) + ')');
                const word f("localForces" + jn + "[pressure,viscous,porous]");
                const word m("localMoments" + jn + "[pressure,viscous,porous]");

                file(i)<< tab << f << tab << m;
            }
        }
    }
    else
    {
        FatalErrorIn("void Foam::forces::writeFileHeader(const label)")
            << "Unhandled file index: " << i
            << abort(FatalError);
    }

    file(i)<< endl;
}


void Foam::forces::initialise()
{
    if (initialised_ || !active_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!obr_.foundObject<volVectorField>(fDName_))
        {
            active_ = false;
            WarningIn("void Foam::forces::initialise()")
                << "Could not find " << fDName_ << " in database." << nl
                << "    De-activating forces."
                << endl;
        }
    }
    else
    {
        if
        (
            !obr_.foundObject<volVectorField>(UName_)
         || !obr_.foundObject<volScalarField>(pName_)
         || (
                rhoName_ != "rhoInf"
             && !obr_.foundObject<volScalarField>(rhoName_)
            )
        )
        {
            active_ = false;

            WarningIn("void Foam::forces::initialise()")
                << "Could not find " << UName_ << ", " << pName_;

            if (rhoName_ != "rhoInf")
            {
                Info<< " or " << rhoName_;
            }

            Info<< " in database." << nl
                << "    De-activating forces." << endl;
        }
    }

    initialised_ = true;
}


Foam::tmp<Foam::volSymmTensorField> Foam::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (obr_.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb =
            obr_.lookupObject<cmpTurbModel>(cmpTurbModel::typeName);

        return turb.devRhoReff();
    }
    else if (obr_.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const incompressible::turbulenceModel& turb =
            obr_.lookupObject<icoTurbModel>(icoTurbModel::typeName);

        return rho()*turb.devReff();
    }
    else if (obr_.foundObject<fluidThermo>(fluidThermo::typeName))
    {
        const fluidThermo& thermo =
            obr_.lookupObject<fluidThermo>(fluidThermo::typeName);

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorIn("forces::devRhoReff()")
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::forces::mu() const
{
    if (obr_.foundObject<fluidThermo>("thermophysicalProperties"))
    {
        const fluidThermo& thermo =
             obr_.lookupObject<fluidThermo>("thermophysicalProperties");

        return thermo.mu();
    }
    else if
    (
        obr_.foundObject<transportModel>("transportProperties")
    )
    {
        const transportModel& laminarT =
            obr_.lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (obr_.foundObject<dictionary>("transportProperties"))
    {
        const dictionary& transportProperties =
             obr_.lookupObject<dictionary>("transportProperties");

        dimensionedScalar nu(transportProperties.lookup("nu"));

        return rho()*nu;
    }
    else
    {
        FatalErrorIn("forces::mu()")
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


Foam::scalar Foam::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorIn("forces::rho(const volScalarField& p)")
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


void Foam::forces::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fT);
        moment_[2][0] += sum(Md^fP);
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label binI = floor(dd[i]/binDx_);
            force_[0][binI] += fN[i];
            force_[1][binI] += fT[i];
            force_[2][binI] += fP[i];
            moment_[0][binI] += Md[i]^fN[i];
            moment_[1][binI] += Md[i]^fT[i];
            moment_[2][binI] += Md[i]^fP[i];
        }
    }
}


void Foam::forces::writeForces()
{
    Info(log_)
        << type() << " " << name_ << " output:" << nl
        << "    sum of forces:" << nl
        << "        pressure : " << sum(force_[0]) << nl
        << "        viscous  : " << sum(force_[1]) << nl
        << "        porous   : " << sum(force_[2]) << nl
        << "    sum of moments:" << nl
        << "        pressure : " << sum(moment_[0]) << nl
        << "        viscous  : " << sum(moment_[1]) << nl
        << "        porous   : " << sum(moment_[2])
        << endl;

    file(0) << obr_.time().value() << tab << setw(1) << '('
        << sum(force_[0]) << setw(1) << ' '
        << sum(force_[1]) << setw(1) << ' '
        << sum(force_[2]) << setw(3) << ") ("
        << sum(moment_[0]) << setw(1) << ' '
        << sum(moment_[1]) << setw(1) << ' '
        << sum(moment_[2]) << setw(1) << ')'
        << endl;

    if (localSystem_)
    {
        vectorField localForceN(coordSys_.localVector(force_[0]));
        vectorField localForceT(coordSys_.localVector(force_[1]));
        vectorField localForceP(coordSys_.localVector(force_[2]));
        vectorField localMomentN(coordSys_.localVector(moment_[0]));
        vectorField localMomentT(coordSys_.localVector(moment_[1]));
        vectorField localMomentP(coordSys_.localVector(moment_[2]));

        file(0) << obr_.time().value() << tab << setw(1) << '('
            << sum(localForceN) << setw(1) << ' '
            << sum(localForceT) << setw(1) << ' '
            << sum(localForceP) << setw(3) << ") ("
            << sum(localMomentN) << setw(1) << ' '
            << sum(localMomentT) << setw(1) << ' '
            << sum(localMomentP) << setw(1) << ')'
            << endl;
    }
}


void Foam::forces::writeBins()
{
    if (nBin_ == 1)
    {
        return;
    }

    List<Field<vector> > f(force_);
    List<Field<vector> > m(moment_);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];

            m[0][i] += m[0][i-1];
            m[1][i] += m[1][i-1];
            m[2][i] += m[2][i-1];
        }
    }

    file(1) << obr_.time().value();

    forAll(f[0], i)
    {
        file(1)
            << tab << setw(1) << '('
            << f[0][i] << setw(1) << ' '
            << f[1][i] << setw(1) << ' '
            << f[2][i] << setw(3) << ") ("
            << m[0][i] << setw(1) << ' '
            << m[1][i] << setw(1) << ' '
            << m[2][i] << setw(1) << ')';
    }

    if (localSystem_)
    {
        List<Field<vector> > lf(3);
        List<Field<vector> > lm(3);
        lf[0] = coordSys_.localVector(force_[0]);
        lf[1] = coordSys_.localVector(force_[1]);
        lf[2] = coordSys_.localVector(force_[2]);
        lm[0] = coordSys_.localVector(moment_[0]);
        lm[1] = coordSys_.localVector(moment_[1]);
        lm[2] = coordSys_.localVector(moment_[2]);

        if (binCumulative_)
        {
            for (label i = 1; i < lf[0].size(); i++)
            {
                lf[0][i] += lf[0][i-1];
                lf[1][i] += lf[1][i-1];
                lf[2][i] += lf[2][i-1];
                lm[0][i] += lm[0][i-1];
                lm[1][i] += lm[1][i-1];
                lm[2][i] += lm[2][i-1];
            }
        }

        forAll(lf[0], i)
        {
            file(1)
                << tab << setw(1) << '('
                << lf[0][i] << setw(1) << ' '
                << lf[1][i] << setw(1) << ' '
                << lf[2][i] << setw(3) << ") ("
                << lm[0][i] << setw(1) << ' '
                << lm[1][i] << setw(1) << ' '
                << lm[2][i] << setw(1) << ')';
        }
    }

    file(1) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    functionObjectFile(obr, name, createFileNames(dict)),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    force_(3),
    moment_(3),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(VGREAT),
    pRef_(0),
    coordSys_(),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(vector::zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        if (readFields)
        {
            read(dict);
            Info<< endl;
        }
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "Foam::forces::forces"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

}


Foam::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const labelHashSet& patchSet,
    const word& pName,
    const word& UName,
    const word& rhoName,
    const scalar rhoInf,
    const scalar pRef,
    const coordinateSystem& coordSys
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    force_(3),
    moment_(3),
    patchSet_(patchSet),
    pName_(pName),
    UName_(UName),
    rhoName_(rhoName),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(rhoInf),
    pRef_(pRef),
    coordSys_(coordSys),
    localSystem_(false),
    porosity_(false),
    nBin_(1),
    binDir_(vector::zero),
    binDx_(0.0),
    binMin_(GREAT),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    forAll(force_, i)
    {
        force_[i].setSize(nBin_);
        moment_[i].setSize(nBin_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::forces::~forces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::forces::read(const dictionary& dict)
{
    if (active_)
    {
        initialised_ = false;

        log_ = dict.lookupOrDefault<Switch>("log", false);

        Info(log_)<< type() << " " << name_ << ":" << nl;

        directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

        if (directForceDensity_)
        {
            // Optional entry for fDName
            fDName_ = dict.lookupOrDefault<word>("fDName", "fD");
        }
        else
        {
            // Optional entries U and p
            pName_ = dict.lookupOrDefault<word>("pName", "p");
            UName_ = dict.lookupOrDefault<word>("UName", "U");
            rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");

            // Reference density needed for incompressible calculations
            rhoRef_ = readScalar(dict.lookup("rhoInf"));

            // Reference pressure, 0 by default
            pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
        }

        coordSys_.clear();

        // Centre of rotation for moment calculations
        // specified directly, from coordinate system, or implicitly (0 0 0)
        if (!dict.readIfPresent<point>("CofR", coordSys_.origin()))
        {
            coordSys_ = coordinateSystem(obr_, dict);
            localSystem_ = true;
        }

        dict.readIfPresent("porosity", porosity_);
        if (porosity_)
        {
            Info(log_)<< "    Including porosity effects" << endl;
        }
        else
        {
            Info(log_)<< "    Not including porosity effects" << endl;
        }

        if (dict.found("binData"))
        {
            const dictionary& binDict(dict.subDict("binData"));
            binDict.lookup("nBin") >> nBin_;

            if (nBin_ < 0)
            {
                FatalIOErrorIn
                (
                    "void Foam::forces::read(const dictionary&)", dict
                )   << "Number of bins (nBin) must be zero or greater"
                    << exit(FatalIOError);
            }
            else if ((nBin_ == 0) || (nBin_ == 1))
            {
                nBin_ = 1;
                forAll(force_, i)
                {
                    force_[i].setSize(1);
                    moment_[i].setSize(1);
                }
            }

            if (nBin_ > 1)
            {
                binDict.lookup("direction") >> binDir_;
                binDir_ /= mag(binDir_);

                binMin_ = GREAT;
                scalar binMax = -GREAT;
                forAllConstIter(labelHashSet, patchSet_, iter)
                {
                    label patchI = iter.key();
                    const polyPatch& pp = pbm[patchI];
                    scalarField d(pp.faceCentres() & binDir_);
                    binMin_ = min(min(d), binMin_);
                    binMax = max(max(d), binMax);
                }
                reduce(binMin_, minOp<scalar>());
                reduce(binMax, maxOp<scalar>());

                // slightly boost binMax so that region of interest is fully
                // within bounds
                binMax = 1.0001*(binMax - binMin_) + binMin_;

                binDx_ = (binMax - binMin_)/scalar(nBin_);

                // create the bin points used for writing
                binPoints_.setSize(nBin_);
                forAll(binPoints_, i)
                {
                    binPoints_[i] = (i + 0.5)*binDir_*binDx_;
                }

                binDict.lookup("cumulative") >> binCumulative_;

                // allocate storage for forces and moments
                forAll(force_, i)
                {
                    force_[i].setSize(nBin_);
                    moment_[i].setSize(nBin_);
                }
            }
        }

        if (nBin_ == 1)
        {
            // allocate storage for forces and moments
            force_[0].setSize(1);
            force_[1].setSize(1);
            force_[2].setSize(1);
            moment_[0].setSize(1);
            moment_[1].setSize(1);
            moment_[2].setSize(1);
        }
    }
}


void Foam::forces::execute()
{
    // Do nothing - only valid on write
}


void Foam::forces::end()
{
    // Do nothing - only valid on write
}


void Foam::forces::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::forces::write()
{
    calcForcesMoment();

    if (!active_)
    {
        return;
    }

    if (Pstream::master())
    {
        functionObjectFile::write();

        writeForces();

        writeBins();

        Info(log_)<< endl;
    }
}


void Foam::forces::calcForcesMoment()
{
    initialise();

    if (!active_)
    {
        return;
    }

    force_[0] = vector::zero;
    force_[1] = vector::zero;
    force_[2] = vector::zero;

    moment_[0] = vector::zero;
    moment_[1] = vector::zero;
    moment_[2] = vector::zero;

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const fvMesh& mesh = fD.mesh();

        const surfaceVectorField::GeometricBoundaryField& Sfb =
            mesh.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            vectorField Md
            (
                mesh.C().boundaryField()[patchI] - coordSys_.origin()
            );

            scalarField sA(mag(Sfb[patchI]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN
            (
                Sfb[patchI]/sA
               *(
                    Sfb[patchI] & fD.boundaryField()[patchI]
                )
            );

            // Tangential force (total force minus normal fN)
            vectorField fT(sA*fD.boundaryField()[patchI] - fN);

            //- Porous force
            vectorField fP(Md.size(), vector::zero);

            applyBins(Md, fN, fT, fP, mesh.C().boundaryField()[patchI]);
        }
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const fvMesh& mesh = U.mesh();

        const surfaceVectorField::GeometricBoundaryField& Sfb =
            mesh.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::GeometricBoundaryField& devRhoReffb
            = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            vectorField Md
            (
                mesh.C().boundaryField()[patchI] - coordSys_.origin()
            );

            vectorField fN
            (
                rho(p)*Sfb[patchI]*(p.boundaryField()[patchI] - pRef)
            );

            vectorField fT(Sfb[patchI] & devRhoReffb[patchI]);

            vectorField fP(Md.size(), vector::zero);

            applyBins(Md, fN, fT, fP, mesh.C().boundaryField()[patchI]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const fvMesh& mesh = U.mesh();

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningIn("void Foam::forces::calcForcesMoment()")
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            // non-const access required if mesh is changing
            porosityModel& pm = const_cast<porosityModel&>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zoneI = cellZoneIDs[i];
                const cellZone& cZone = mesh.cellZones()[zoneI];

                const vectorField d(mesh.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - coordSys_.origin());

                const vectorField fDummy(Md.size(), vector::zero);

                applyBins(Md, fDummy, fDummy, fP, d);
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::forces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


Foam::vector Foam::forces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]);
}


// ************************************************************************* //
