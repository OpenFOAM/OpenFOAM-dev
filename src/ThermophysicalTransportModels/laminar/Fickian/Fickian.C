/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

    mixtureDiffusionCoefficients_(true),

    DFuncs_(this->thermo().composition().species().size()),

    DmFuncs_(this->thermo().composition().species().size()),

    DTFuncs_
    (
        this->coeffDict_.found("DT")
      ? this->thermo().composition().species().size()
      : 0
    ),

    Dm_(this->thermo().composition().species().size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Fickian<BasicThermophysicalTransportModel>::read()
{
    if
    (
        BasicThermophysicalTransportModel::read()
    )
    {
        const basicSpecieMixture& composition = this->thermo().composition();
        const speciesTable& species = composition.species();

        this->coeffDict_.lookup("mixtureDiffusionCoefficients")
            >> mixtureDiffusionCoefficients_;

        if (mixtureDiffusionCoefficients_)
        {
            const dictionary& Ddict = this->coeffDict_.subDict("Dm");

            forAll(species, i)
            {
                DmFuncs_.set
                (
                    i,
                    Function2<scalar>::New(species[i], Ddict).ptr()
                );
            }
        }
        else
        {
            const dictionary& Ddict = this->coeffDict_.subDict("D");

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
                            Function2<scalar>::New(Dname, Ddict).ptr()
                        );
                    }
                }
            }
        }

        // Optionally read the List of specie Soret thermal diffusion
        // coefficient functions
        if (this->coeffDict_.found("DT"))
        {
            const dictionary& DTdict = this->coeffDict_.subDict("DT");

            forAll(species, i)
            {
                DTFuncs_.set
                (
                    i,
                    Function2<scalar>::New(species[i], DTdict).ptr()
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
    const basicSpecieMixture& composition = this->thermo().composition();

    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *Dm_[composition.index(Yi)]
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> Fickian<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();

    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dm_[composition.index(Yi)].boundaryField()[patchi];
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

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

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
            if (i != composition.defaultSpecie())
            {
                const volScalarField hi
                (
                    composition.Hs(i, this->thermo().p(), this->thermo().T())
                );

                const surfaceScalarField ji(this->j(Y[i]));
                sumJ += ji;

                sumJh += ji*fvc::interpolate(hi);
            }
        }

        {
            const label i = composition.defaultSpecie();

            const volScalarField hi
            (
                composition.Hs(i, this->thermo().p(), this->thermo().T())
            );

            sumJh -= sumJ*fvc::interpolate(hi);
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

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();

    tmpDivq.ref() -=
        correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));

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
        if (i != composition.defaultSpecie())
        {
            const volScalarField hi
            (
                composition.Hs(i, this->thermo().p(), this->thermo().T())
            );

            const surfaceScalarField ji(this->j(Y[i]));
            sumJ += ji;

            sumJh += ji*fvc::interpolate(hi);
        }
    }

    {
        const label i = composition.defaultSpecie();

        const volScalarField hi
        (
            composition.Hs(i, this->thermo().p(), this->thermo().T())
        );

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
        const basicSpecieMixture& composition = this->thermo().composition();
        const volScalarField& p = this->thermo().p();
        const volScalarField& T = this->thermo().T();

        return
            BasicThermophysicalTransportModel::j(Yi)
          - fvc::interpolate
            (
                evaluate
                (
                    DTFuncs_[composition.index(Yi)],
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
tmp<fvScalarMatrix> Fickian<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    if (DTFuncs_.size())
    {
        const basicSpecieMixture& composition = this->thermo().composition();
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
                        DTFuncs_[composition.index(Yi)],
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
void Fickian<BasicThermophysicalTransportModel>::correct()
{
    BasicThermophysicalTransportModel::correct();

    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    if (mixtureDiffusionCoefficients_)
    {
        forAll(Y, i)
        {
            Dm_.set(i, evaluate(DmFuncs_[i], dimViscosity, p, T));
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
                dimless/dimViscosity/Wm.dimensions()
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
                           dimensionedScalar
                           (
                               "Wj",
                               Wm.dimensions(),
                               composition.Wi(j)
                           )
                          *(
                               i < j
                             ? evaluate(DFuncs_[i][j], dimViscosity, p, T)
                             : evaluate(DFuncs_[j][i], dimViscosity, p, T)
                           )
                       );
                }
            }

            Dm_.set
            (
                i,
                (
                    1/Wm
                  - Y[i]
                   /dimensionedScalar("Wi", Wm.dimensions(), composition.Wi(i))
                )/max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
