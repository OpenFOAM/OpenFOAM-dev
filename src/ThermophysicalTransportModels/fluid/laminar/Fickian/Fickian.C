/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "Fickian.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "Function2Evaluate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
void Fickian<BasicThermophysicalTransportModel>::updateDm() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    Dm_.setSize(Y.size());

    if (mixtureDiffusionCoefficients_)
    {
        forAll(Y, i)
        {
            Dm_.set(i, evaluate(DmFuncs_[i], dimKinematicViscosity, p, T));
        }
    }
    else
    {
        const volScalarField Wm(this->thermo().W());
        volScalarField sumXbyD
        (
            volScalarField::New
            (
                "sumXbyD",
                T.mesh(),
                dimless/dimKinematicViscosity/Wm.dimensions()
            )
        );

        forAll(Dm_, i)
        {
            sumXbyD = Zero;

            forAll(Y, j)
            {
                if (j != i)
                {
                    sumXbyD +=
                        Y[j]
                       /(
                            this->thermo().Wi(j)
                           *(
                                i < j
                              ? evaluate
                                (
                                    DFuncs_[i][j],
                                    dimKinematicViscosity,
                                    p,
                                    T
                                )
                              : evaluate
                                (
                                    DFuncs_[j][i],
                                    dimKinematicViscosity,
                                    p,
                                    T
                                )
                           )
                       );
                }
            }

            Dm_.set
            (
                i,
                (1/Wm - Y[i]/this->thermo().Wi(i))
               /max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small))
            );
        }
    }
}


template<class BasicThermophysicalTransportModel>
const PtrList<volScalarField>&
Fickian<BasicThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }

    return Dm_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Fickian<BasicThermophysicalTransportModel>::Fickian
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel
    (
        type,
        momentumTransport,
        thermo
    ),

    TopoChangeableMeshObject(*this),

    mixtureDiffusionCoefficients_(true),

    DFuncs_(this->thermo().species().size()),

    DmFuncs_(this->thermo().species().size()),

    DTFuncs_
    (
        this->coeffDict().found("DT")
      ? this->thermo().species().size()
      : 0
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Fickian<BasicThermophysicalTransportModel>::read()
{
    if (BasicThermophysicalTransportModel::read())
    {
        const speciesTable& species = this->thermo().species();

        const dictionary& coeffDict = this->coeffDict();

        coeffDict.lookup("mixtureDiffusionCoefficients")
            >> mixtureDiffusionCoefficients_;

        if (mixtureDiffusionCoefficients_)
        {
            const dictionary& Ddict = coeffDict.subDict("Dm");

            forAll(species, i)
            {
                DmFuncs_.set
                (
                    i,
                    Function2<scalar>::New
                    (
                        species[i],
                        dimPressure,
                        dimTemperature,
                        dimKinematicViscosity,
                        Ddict
                    ).ptr()
                );
            }
        }
        else
        {
            const dictionary& Ddict = coeffDict.subDict("D");

            // Read the array of specie binary mass diffusion coefficient
            // functions
            forAll(species, i)
            {
                DFuncs_[i].setSize(species.size());

                forAll(species, j)
                {
                    if (j >= i)
                    {
                        const word nameij(species[i] + '-' + species[j]);
                        const word nameji(species[j] + '-' + species[i]);

                        word Dname;

                        if (Ddict.found(nameij) && Ddict.found(nameji))
                        {
                            if (i != j)
                            {
                                WarningInFunction
                                    << "Binary mass diffusion coefficients "
                                       "for Both " << nameij
                                    << " and " << nameji << " provided, using "
                                    << nameij << endl;
                            }

                            Dname = nameij;
                        }
                        else if (Ddict.found(nameij))
                        {
                            Dname = nameij;
                        }
                        else if (Ddict.found(nameji))
                        {
                            Dname = nameji;
                        }
                        else
                        {
                            FatalIOErrorInFunction(Ddict)
                                << "Binary mass diffusion coefficient for pair "
                                << nameij << " or " << nameji << " not provided"
                                << exit(FatalIOError);
                        }

                        DFuncs_[i].set
                        (
                            j,
                            Function2<scalar>::New
                            (
                                Dname,
                                dimPressure,
                                dimTemperature,
                                dimKinematicViscosity,
                                Ddict
                            ).ptr()
                        );
                    }
                }
            }
        }

        // Optionally read the List of specie Soret thermal diffusion
        // coefficient functions
        if (coeffDict.found("DT"))
        {
            const dictionary& DTdict = coeffDict.subDict("DT");

            forAll(species, i)
            {
                DTFuncs_.set
                (
                    i,
                    Function2<scalar>::New
                    (
                        species[i],
                        dimPressure,
                        dimTemperature,
                        dimDynamicViscosity,
                        DTdict
                    ).ptr()
                );
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicThermophysicalTransportModel>
tmp<volScalarField> Fickian<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *Dm()[this->thermo().specieIndex(Yi)]
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fickian<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dm()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fickian<BasicThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    if (Y.size())
    {
        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(sumJ.dimensions()*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const volScalarField hi(this->thermo().hsi(i, p, T));

                const surfaceScalarField ji(this->j(Y[i]));

                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const volScalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*fvc::interpolate(hi);
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fickian<BasicThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    tmp<scalarField> tmpq
    (
      - (
            this->alpha().boundaryField()[patchi]
           *this->kappaEff(patchi)
           *this->thermo().T().boundaryField()[patchi].snGrad()
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const scalarField& p = this->thermo().p().boundaryField()[patchi];
    const scalarField& T = this->thermo().T().boundaryField()[patchi];

    if (Y.size())
    {
        scalarField sumJ(tmpq->size(), scalar(0));
        scalarField sumJh(tmpq->size(), scalar(0));

        forAll(Y, i)
        {
            if (i != this->thermo().defaultSpecie())
            {
                const scalarField hi(this->thermo().hsi(i, p, T));

                const scalarField ji(this->j(Y[i], patchi));

                sumJ += ji;

                sumJh += ji*hi;
            }
        }

        {
            const label i = this->thermo().defaultSpecie();

            const scalarField hi(this->thermo().hsi(i, p, T));

            sumJh -= sumJ*hi;
        }

        tmpq.ref() += sumJh;
    }

    return tmpq;
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> Fickian<BasicThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    tmp<fvScalarMatrix> tmpDivq
    (
        fvm::Su
        (
            -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
            he
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    tmpDivq.ref() -=
        fvm::laplacianCorrection(this->alpha()*this->alphaEff(), he);

    surfaceScalarField sumJ
    (
        surfaceScalarField::New
        (
            "sumJ",
            he.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        )
    );

    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "sumJh",
            he.mesh(),
            dimensionedScalar(sumJ.dimensions()*he.dimensions(), 0)
        )
    );

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));

            const surfaceScalarField ji(this->j(Y[i]));

            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = this->thermo().defaultSpecie();

        const volScalarField hi(this->thermo().hsi(i, p, T));

        sumJh -= sumJ*fvc::interpolate(hi);
    }

    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());

    return tmpDivq;
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField> Fickian<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    if (DTFuncs_.size())
    {
        const volScalarField& p = this->thermo().p();
        const volScalarField& T = this->thermo().T();

        return
            BasicThermophysicalTransportModel::j(Yi)
          - fvc::interpolate
            (
                evaluate
                (
                    DTFuncs_[this->thermo().specieIndex(Yi)],
                    dimDynamicViscosity,
                    p,
                    T
                )
            )
           *fvc::snGrad(T)/fvc::interpolate(T);
    }
    else
    {
        return BasicThermophysicalTransportModel::j(Yi);
    }
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fickian<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    if (DTFuncs_.size())
    {
        const scalarField& p = this->thermo().p().boundaryField()[patchi];
        const scalarField& T = this->thermo().T().boundaryField()[patchi];

        return
            BasicThermophysicalTransportModel::j(Yi, patchi)
          - DTFuncs_[this->thermo().specieIndex(Yi)].value(p, T)
           *this->thermo().T().boundaryField()[patchi].snGrad()/T;
    }
    else
    {
        return BasicThermophysicalTransportModel::j(Yi, patchi);
    }
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> Fickian<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    if (DTFuncs_.size())
    {
        const volScalarField& p = this->thermo().p();
        const volScalarField& T = this->thermo().T();

        return
            BasicThermophysicalTransportModel::divj(Yi)
          - fvc::div
            (
                fvc::interpolate
                (
                    evaluate
                    (
                        DTFuncs_[this->thermo().specieIndex(Yi)],
                        dimDynamicViscosity,
                        p,
                        T
                    )
                )
               *fvc::snGrad(T)/fvc::interpolate(T)
               *T.mesh().magSf()
            );
    }
    else
    {
        return BasicThermophysicalTransportModel::divj(Yi);
    }
}


template<class BasicThermophysicalTransportModel>
void Fickian<BasicThermophysicalTransportModel>::predict()
{
    BasicThermophysicalTransportModel::predict();
    updateDm();
}


template<class BasicThermophysicalTransportModel>
bool Fickian<BasicThermophysicalTransportModel>::movePoints()
{
    return true;
}


template<class BasicThermophysicalTransportModel>
void Fickian<BasicThermophysicalTransportModel>::topoChange
(
    const polyTopoChangeMap& map
)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}


template<class BasicThermophysicalTransportModel>
void Fickian<BasicThermophysicalTransportModel>::mapMesh
(
    const polyMeshMap& map
)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}


template<class BasicThermophysicalTransportModel>
void Fickian<BasicThermophysicalTransportModel>::distribute
(
    const polyDistributionMap& map
)
{
    // Delete the cached Dm, will be re-created in predict
    Dm_.clear();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
