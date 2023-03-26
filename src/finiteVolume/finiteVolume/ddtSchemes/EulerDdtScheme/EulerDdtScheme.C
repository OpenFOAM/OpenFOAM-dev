/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "EulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+dt.name()+')');

    if (mesh().moving())
    {
        tmp<VolField<Type>> tdtdt
        (
            VolField<Type>::New
            (
                ddtName,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    Zero
                )
            )
        );

        tdtdt.ref().primitiveFieldRef() =
            rDeltaT.value()*dt.value()*(1.0 - mesh().Vsc0()/mesh().Vsc());

        return tdtdt;
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                Zero
            ),
            calculatedFvPatchField<Type>::typeName
        );
    }
}


template<class Type>
tmp<VolField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+vf.name()+')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*
            (
                vf()
              - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()*
            (
                vf.boundaryField() - vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*(vf - vf.oldTime())
        );
    }
}


template<class Type>
tmp<VolField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*rho*
            (
                vf()
              - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()*rho.value()*
            (
                vf.boundaryField() - vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*rho*(vf - vf.oldTime())
        );
    }
}


template<class Type>
tmp<VolField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*
            (
                rho()*vf()
              - rho.oldTime()()
               *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()*
            (
                rho.boundaryField()*vf.boundaryField()
              - rho.oldTime().boundaryField()
               *vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
        );
    }
}


template<class Type>
tmp<VolField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+alpha.name()+','+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT*
            (
                alpha()
               *rho()
               *vf()

              - alpha.oldTime()()
               *rho.oldTime()()
               *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
            ),
            rDeltaT.value()*
            (
                alpha.boundaryField()
               *rho.boundaryField()
               *vf.boundaryField()

              - alpha.oldTime().boundaryField()
               *rho.oldTime().boundaryField()
               *vf.oldTime().boundaryField()
            )
        );
    }
    else
    {
        return VolField<Type>::New
        (
            ddtName,
            rDeltaT
           *(
               alpha*rho*vf
             - alpha.oldTime()*rho.oldTime()*vf.oldTime()
            )
        );
    }
}


template<class Type>
tmp<SurfaceField<Type>>
EulerDdtScheme<Type>::fvcDdt
(
    const SurfaceField<Type>& sf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+sf.name()+')');

    return SurfaceField<Type>::New
    (
        ddtName,
        rDeltaT*(sf - sf.oldTime())
    );
}


template<class Type>
tmp<fvMatrix<Type>>
EulerDdtScheme<Type>::fvmDdt
(
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
EulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.value()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
EulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
EulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    const scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() =
        rDeltaT*alpha.primitiveField()*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtUfCorr
(
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
    (
        phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr) *rDeltaT*phiCorr
    );
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    fluxFieldType phiCorr
    (
        phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
       *rDeltaT*phiCorr
    );
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& rhoUf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if
    (
        U.dimensions() == dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        VolField<Type> rhoU0(rho.oldTime()*U.oldTime());
        fluxFieldType phiUf0(mesh().Sf() & rhoUf.oldTime());
        fluxFieldType phiCorr(phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0));

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + rhoUf.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr, rho.oldTime())
           *rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        fluxFieldType phiUf0(mesh().Sf() & rhoUf.oldTime());
        fluxFieldType phiCorr
        (
            phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + rhoUf.name() + ')',
            this->fvcDdtPhiCoeff
            (
                U.oldTime(),
                phiUf0,
                phiCorr,
                rho.oldTime()
            )*rDeltaT*phiCorr
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of rhoUf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimFlux
    )
    {
        VolField<Type> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), rhoU0)
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
            this->fvcDdtPhiCoeff
            (
                rhoU0,
                phi.oldTime(),
                phiCorr,
                rho.oldTime()
            )*rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimFlux
    )
    {
        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
            this->fvcDdtPhiCoeff
            (
                U.oldTime(),
                phi.oldTime(),
                phiCorr,
                rho.oldTime()
            )*rDeltaT*phiCorr
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if (U.dimensions() == dimVelocity && Uf.dimensions() == dimVelocity)
    {
        const volScalarField alphaRho0(alpha.oldTime()*rho.oldTime());

        return fluxFieldType::New
        (
            "ddtCorr("
          + alpha.name() + rho.name() + ',' + U.name() + ',' + Uf.name()
          + ')',
            this->fvcDdtPhiCoeff(U.oldTime(), mesh().Sf() & Uf.oldTime())
           *rDeltaT
           *(
                mesh().Sf()
              & (
                  fvc::interpolate(alphaRho0)*Uf.oldTime()
                - fvc::interpolate(alphaRho0*U.oldTime())
                )
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename EulerDdtScheme<Type>::fluxFieldType>
EulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if (U.dimensions() == dimVelocity && phi.dimensions() == dimFlux)
    {
        const volScalarField alphaRho0(alpha.oldTime()*rho.oldTime());

        return fluxFieldType::New
        (
            "ddtCorr("
          + alpha.name() + rho.name() + ',' + U.name() + ',' + phi.name()
          + ')',
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *rDeltaT
           *(
                fvc::interpolate(alphaRho0)*phi.oldTime()
               -fvc::dotInterpolate(mesh().Sf(), alphaRho0*U.oldTime())
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> EulerDdtScheme<Type>::meshPhi
(
    const VolField<Type>&
)
{
    return mesh().phi();
}


template<class Type>
tmp<scalarField> EulerDdtScheme<Type>::meshPhi
(
    const VolField<Type>&,
    const label patchi
)
{
    return mesh().phi().boundaryField()[patchi];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
