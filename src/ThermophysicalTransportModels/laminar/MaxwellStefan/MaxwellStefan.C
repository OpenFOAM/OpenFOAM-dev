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

#include "MaxwellStefan.H"
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
MaxwellStefan<BasicThermophysicalTransportModel>::MaxwellStefan
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

    DFuncs_(this->thermo().composition().species().size()),

    DTFuncs_
    (
        this->coeffDict_.found("DT")
      ? this->thermo().composition().species().size()
      : 0
    ),

    Dii_(this->thermo().composition().species().size()),
    jexp_(this->thermo().composition().species().size()),

    W(this->thermo().composition().species().size()),

    YPtrs(W.size()),
    DijPtrs(W.size()),

    Y(W.size()),
    X(W.size()),
    DD(W.size()),
    A(W.size() - 1),
    B(A.m()),
    invA(A.m()),
    D(W.size())
{
    const basicSpecieMixture& composition = this->thermo().composition();

    // Set the molecular weights of the species
    forAll(W, i)
    {
        W[i] = composition.Wi(i);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool MaxwellStefan<BasicThermophysicalTransportModel>::read()
{
    if
    (
        BasicThermophysicalTransportModel::read()
    )
    {
        const basicSpecieMixture& composition = this->thermo().composition();
        const speciesTable& species = composition.species();

        const dictionary& Ddict = this->coeffDict_.subDict("D");

        // Read the array of specie binary mass diffusion coefficient functions
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
                                   "for both " << nameij << " and " << nameji
                                << " provided, using " << nameij << endl;
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
                            << "Binary mass diffusion coefficients for pair "
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
tmp<volScalarField> MaxwellStefan<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();

    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()*Dii_[composition.index(Yi)]
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField> MaxwellStefan<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();

    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *Dii_[composition.index(Yi)].boundaryField()[patchi];
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
MaxwellStefan<BasicThermophysicalTransportModel>::q() const
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
    const label d = composition.defaultSpecie();

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
            if (i != d)
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
            const label i = d;

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
tmp<fvScalarMatrix> MaxwellStefan<BasicThermophysicalTransportModel>::divq
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
    const label d = composition.defaultSpecie();

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
        if (i != d)
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
        const label i = d;

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
tmp<surfaceScalarField> MaxwellStefan<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();
    return
        BasicThermophysicalTransportModel::j(Yi)
      + jexp_[composition.index(Yi)];
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix> MaxwellStefan<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const basicSpecieMixture& composition = this->thermo().composition();
    return
        BasicThermophysicalTransportModel::divj(Yi)
      + fvc::div(jexp_[composition.index(Yi)]*Yi.mesh().magSf());
}


template<class BasicThermophysicalTransportModel>
void MaxwellStefan<BasicThermophysicalTransportModel>::
transformDiffusionCoefficient()
{
    const basicSpecieMixture& composition = this->thermo().composition();
    const label d = composition.defaultSpecie();

    // Calculate the molecular weight of the mixture and the mole fractions
    scalar Wm = 0;

    forAll(W, i)
    {
        X[i] = Y[i]/W[i];
        Wm += X[i];
    }

    Wm = 1/Wm;
    X *= Wm;

    // i counter for the specie sub-system without the default specie
    label is = 0;

    // Calculate the A and B matrices from the binary mass diffusion
    // coefficients and specie mole fractions
    forAll(X, i)
    {
        if (i != d)
        {
            A(is, is) = -X[i]*Wm/(DD(i, d)*W[d]);
            B(is, is) = -(X[i]*Wm/W[d] + (1 - X[i])*Wm/W[i]);

            // j counter for the specie sub-system without the default specie
            label js = 0;

            forAll(X, j)
            {
                if (j != i)
                {
                    A(is, is) -= X[j]*Wm/(DD(i, j)*W[i]);

                    if (j != d)
                    {
                        A(is, js) =
                            X[i]*(Wm/(DD(i, j)*W[j]) - Wm/(DD(i, d)*W[d]));

                        B(is, js) = X[i]*(Wm/W[j] - Wm/W[d]);
                    }
                }

                if (j != d)
                {
                    js++;
                }
            }

            is++;
        }
    }

    // LU decompose A and invert
    A.decompose();
    A.inv(invA);

    // Calculate the generalised Fick's law diffusion coefficients
    multiply(D, invA, B);
}


template<class BasicThermophysicalTransportModel>
void MaxwellStefan<BasicThermophysicalTransportModel>::
transformDiffusionCoefficientFields()
{
    const basicSpecieMixture& composition = this->thermo().composition();
    const label d = composition.defaultSpecie();

    // For each cell or patch face
    forAll(*(YPtrs[0]), pi)
    {
        forAll(W, i)
        {
            // Map YPtrs -> Y
            Y[i] = (*YPtrs[i])[pi];

            // Map DijPtrs -> DD
            forAll(W, j)
            {
                DD(i, j) = (*DijPtrs[i][j])[pi];
            }
        }

        // Transform DD -> D
        transformDiffusionCoefficient();

        // i counter for the specie sub-system without the default specie
        label is = 0;

        forAll(W, i)
        {
            if (i != d)
            {
                // j counter for the specie sub-system
                // without the default specie
                label js = 0;

                // Map D -> DijPtrs
                forAll(W, j)
                {
                    if (j != d)
                    {
                        (*DijPtrs[i][j])[pi] = D(is, js);

                        js++;
                    }
                }

                is++;
            }
        }
    }
}


template<class BasicThermophysicalTransportModel>
void MaxwellStefan<BasicThermophysicalTransportModel>::transform
(
    List<PtrList<volScalarField>>& Dij
)
{
    const basicSpecieMixture& composition = this->thermo().composition();
    const PtrList<volScalarField>& Y = composition.Y();
    const volScalarField& Y0 = Y[0];

    forAll(W, i)
    {
        // Map composition.Y() internal fields -> YPtrs
        YPtrs[i] = &Y[i].primitiveField();

        // Map Dii_ internal fields -> DijPtrs
        DijPtrs[i][i] = &Dii_[i].primitiveFieldRef();

        // Map Dij internal fields -> DijPtrs
        forAll(W, j)
        {
            if (j != i)
            {
                DijPtrs[i][j] = &Dij[i][j].primitiveFieldRef();
            }
        }
    }

    // Transform binary mass diffusion coefficients internal field DijPtrs ->
    // generalised Fick's law diffusion coefficients DijPtrs
    transformDiffusionCoefficientFields();

    forAll(Y0.boundaryField(), patchi)
    {
        forAll(W, i)
        {
            // Map composition.Y() patch fields -> YPtrs
            YPtrs[i] = &Y[i].boundaryField()[patchi];

            // Map Dii_ patch fields -> DijPtrs
            DijPtrs[i][i] = &Dii_[i].boundaryFieldRef()[patchi];

            // Map Dij patch fields -> DijPtrs
            forAll(W, j)
            {
                if (j != i)
                {
                    DijPtrs[i][j] = &Dij[i][j].boundaryFieldRef()[patchi];
                }
            }
        }

        // Transform binary mass diffusion coefficients patch field DijPtrs ->
        // generalised Fick's law diffusion coefficients DijPtrs
        transformDiffusionCoefficientFields();
    }
}


template<class BasicThermophysicalTransportModel>
void MaxwellStefan<BasicThermophysicalTransportModel>::correct()
{
    BasicThermophysicalTransportModel::correct();

    const basicSpecieMixture& composition = this->thermo().composition();
    const label d = composition.defaultSpecie();

    const PtrList<volScalarField>& Y = composition.Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField& rho = this->momentumTransport().rho();

    List<PtrList<volScalarField>> Dij(Y.size());

    // Evaluate the specie binary mass diffusion coefficient functions
    // and initialise the explicit part of the specie mass flux fields
    forAll(Y, i)
    {
        if (i != d)
        {
            if (jexp_.set(i))
            {
                jexp_[i] = Zero;
            }
            else
            {
                jexp_.set
                (
                    i,
                    surfaceScalarField::New
                    (
                        "jexp" + Y[i].name(),
                        Y[i].mesh(),
                        dimensionedScalar(dimensionSet(1, -2, -1, 0, 0), 0)
                    )
                );
            }
        }

        Dii_.set(i, evaluate(DFuncs_[i][i], dimViscosity, p, T));

        Dij[i].setSize(Y.size());

        forAll(Y, j)
        {
            if (j > i)
            {
                Dij[i].set(j, evaluate(DFuncs_[i][j], dimViscosity, p, T));
            }
            else if (j < i)
            {
                Dij[i].set(j, Dij[j][i].clone());
            }
        }
    }

    //- Transform the binary mass diffusion coefficients into the
    //  the generalised Fick's law diffusion coefficients
    transform(Dij);

    // Accumulate the explicit part of the specie mass flux fields
    forAll(Y, j)
    {
        if (j != d)
        {
            const surfaceScalarField snGradYj(fvc::snGrad(Y[j]));

            forAll(Y, i)
            {
                if (i != d && i != j)
                {
                    jexp_[i] -= fvc::interpolate(rho*Dij[i][j])*snGradYj;
                }
            }
        }
    }

    // Optionally add the Soret thermal diffusion contribution to the
    // explicit part of the specie mass flux fields
    if (DTFuncs_.size())
    {
        const surfaceScalarField gradTbyT(fvc::snGrad(T)/fvc::interpolate(T));

        forAll(Y, i)
        {
            if (i != d)
            {
                jexp_[i] -= fvc::interpolate
                (
                    evaluate(DTFuncs_[i], dimDynamicViscosity, p, T)
                )*gradTbyT;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
