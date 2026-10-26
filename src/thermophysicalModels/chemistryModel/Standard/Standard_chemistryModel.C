/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "Standard_chemistryModel.H"
#include "multicomponentMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModels::Standard<ThermoType>::Standard
(
    const fluidMulticomponentThermo& thermo
)
:
    standard
    (
        thermo,
        SpecieThermosRef<ThermoType>
        (
            dynamicCast<const multicomponentMixture<ThermoType>>
            (
                thermo
            ).specieThermos()
        )
    ),
    jacobianType_
    (
        this->found("jacobian")
      ? jacobianTypeNames.read(this->lookup("jacobian"))
      : jacobianType::fast
    ),
    specieThermos_
    (
        dynamicCast<const multicomponentMixture<ThermoType>>
        (
            this->thermo()
        ).specieThermos()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModels::Standard<ThermoType>::~Standard()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryModels::Standard<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& YT,
    const label li,
    scalarField& dYTdt
) const
{
    if (reduced_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YT[i], 0);
        }
    }
    else
    {
        forAll(Y_, i)
        {
            Y_[i] = max(YT[i], 0);
        }
    }

    const scalar T
    (
        maxMin(YT[nSpecie_], reactions_.Tlow(), reactions_.Thigh())
    );

    // Evaluate the mixture density
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]/specieThermos_[i].rho(p_, T);
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        s_.c[i] = rhoM*Y_[i]/specieThermos_[i].W();
    }

    s_.reset(p_, T, li);

    // Evaluate contributions from reactions
    dYTdt = Zero;
    forAll(reactions_, ri)
    {
        if (!reduction_.reactionDisabled(ri))
        {
            reactions_[ri].dNdtByV(s_, dYTdt, reduced_, c2s_);
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTdt[i];
        dYidt *= WiByrhoM;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp
    scalar CpM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        CpM += Y_[i]*specieThermos_[i].Cp(p_, T);
    }

    // dT/dt
    scalar& dTdt = dYTdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dYTdt[i]*specieThermos_[sToc(i)].ha(p_, T);
    }
    dTdt /= CpM;
}


template<class ThermoType>
void Foam::chemistryModels::Standard<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YT,
    const label li,
    scalarField& dYTdt,
    scalarSquareMatrix& J
) const
{
    if (reduced_)
    {
        forAll(sToc_, i)
        {
            Y_[sToc_[i]] = max(YT[i], 0);
        }
    }
    else
    {
        forAll(Y_, i)
        {
            Y_[i] = max(YT[i], 0);
        }
    }

    const scalar T
    (
        maxMin(YT[nSpecie_], reactions_.Tlow(), reactions_.Thigh())
    );

    // Evaluate the specific volumes and mixture density
    scalarField& v = YTWork_.elmt<0>();
    for (label i=0; i<Y_.size(); i++)
    {
        v[i] = 1/specieThermos_[i].rho(p_, T);
    }
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]*v[i];
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        s_.c[i] = rhoM*Y_[i]/specieThermos_[i].W();
    }

    s_.reset(p_, T, li);

    // Evaluate the derivatives of concentration w.r.t. mass fraction
    scalarSquareMatrix& dcdY = YTYTWork_.elmt<0>();
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar rhoMByWi = rhoM/specieThermos_[sToc(i)].W();
        switch (jacobianType_)
        {
            case jacobianType::fast:
                {
                    dcdY(i, i) = rhoMByWi;
                }
                break;
            case jacobianType::exact:
                for (label j=0; j<nSpecie_; j++)
                {
                    dcdY(i, j) =
                        rhoMByWi*((i == j) - rhoM*v[sToc(j)]*Y_[sToc(i)]);
                }
                break;
        }
    }

    // Evaluate the mixture thermal expansion coefficient
    scalar alphavM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        alphavM += Y_[i]*rhoM*v[i]*specieThermos_[i].alphav(p_, T);
    }

    // Evaluate contributions from reactions
    dYTdt = Zero;
    scalarSquareMatrix& ddNdtByVdcT = YTYTWork_.elmt<1>();
    for (label i=0; i<nSpecie_ + 1; i++)
    {
        for (label j=0; j<nSpecie_ + 1; j++)
        {
            ddNdtByVdcT[i][j] = 0;
        }
    }
    forAll(reactions_, ri)
    {
        if (!reduction_.reactionDisabled(ri))
        {
            reactions_[ri].ddNdtByVdcT
            (
                s_,
                dYTdt,
                ddNdtByVdcT,
                reduced_,
                c2s_,
                nSpecie_
            );
        }
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTdt[i];
        dYidt *= WiByrhoM;

        for (label j=0; j<nSpecie_; j++)
        {
            scalar ddNidtByVdYj = 0;
            switch (jacobianType_)
            {
                case jacobianType::fast:
                    {
                        const scalar ddNidtByVdcj = ddNdtByVdcT(i, j);
                        ddNidtByVdYj = ddNidtByVdcj*dcdY(j, j);
                    }
                    break;
                case jacobianType::exact:
                    for (label k=0; k<nSpecie_; k++)
                    {
                        const scalar ddNidtByVdck = ddNdtByVdcT(i, k);
                        ddNidtByVdYj += ddNidtByVdck*dcdY(k, j);
                    }
                    break;
            }

            scalar& ddYidtdYj = J(i, j);
            ddYidtdYj = WiByrhoM*ddNidtByVdYj + rhoM*v[sToc(j)]*dYidt;
        }

        scalar ddNidtByVdT = ddNdtByVdcT(i, nSpecie_);
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddNidtByVdcj = ddNdtByVdcT(i, j);
            ddNidtByVdT -= ddNidtByVdcj*s_.c[sToc(j)]*alphavM;
        }

        scalar& ddYidtdT = J(i, nSpecie_);
        ddYidtdT = WiByrhoM*ddNidtByVdT + alphavM*dYidt;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp and its derivative
    scalarField& Cp = YTWork_.elmt<1>();
    scalar CpM = 0, dCpMdT = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        Cp[i] = specieThermos_[i].Cp(p_, T);
        CpM += Y_[i]*Cp[i];
        dCpMdT += Y_[i]*specieThermos_[i].dCpdT(p_, T);
    }

    // dT/dt
    scalarField& ha = YTWork_.elmt<2>();
    scalar& dTdt = dYTdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        ha[sToc(i)] = specieThermos_[sToc(i)].ha(p_, T);
        dTdt -= dYTdt[i]*ha[sToc(i)];
    }
    dTdt /= CpM;

    // d(dTdt)/dY
    for (label i=0; i<nSpecie_; i++)
    {
        scalar& ddTdtdYi = J(nSpecie_, i);
        ddTdtdYi = 0;
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddYjdtdYi = J(j, i);
            ddTdtdYi -= ddYjdtdYi*ha[sToc(j)];
        }
        ddTdtdYi -= Cp[sToc(i)]*dTdt;
        ddTdtdYi /= CpM;
    }

    // d(dTdt)/dT
    scalar& ddTdtdT = J(nSpecie_, nSpecie_);
    ddTdtdT = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar dYidt = dYTdt[i];
        const scalar ddYidtdT = J(i, nSpecie_);
        ddTdtdT -= dYidt*Cp[sToc(i)] + ddYidtdT*ha[sToc(i)];
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT /= CpM;
}


// ************************************************************************* //
