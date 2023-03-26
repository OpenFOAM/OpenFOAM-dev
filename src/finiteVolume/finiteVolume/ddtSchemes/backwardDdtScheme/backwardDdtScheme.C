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

#include "backwardDdtScheme.H"
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
scalar backwardDdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar backwardDdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
template<class GeoField>
scalar backwardDdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    if (vf.nOldTimes() < 2)
    {
        return great;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
backwardDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+dt.name()+')');

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_();

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

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

        tdtdt.ref().primitiveFieldRef() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().V0() - coefft00*mesh().V00())/mesh().V()
        );

        return tdtdt;
    }
    else
    {
        return tmp<VolField<Type>>
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
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<VolField<Type>>
backwardDdtScheme<Type>::fvcDdt
(
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+vf.name()+')');

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*vf() -
                    (
                        coefft0*vf.oldTime()()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<VolField<Type>>
backwardDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*rho*
                (
                    coefft*vf() -
                    (
                        coefft0*vf.oldTime()()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime()()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<VolField<Type>>
backwardDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*rho()*vf() -
                    (
                        coefft0*rho.oldTime()()
                       *vf.oldTime()()*mesh().V0()
                      - coefft00*rho.oldTime().oldTime()()
                       *vf.oldTime().oldTime()()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*rho.boundaryField()*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<VolField<Type>>
backwardDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const word ddtName("ddt("+alpha.name()+','+rho.name()+','+vf.name()+')');

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*alpha()*rho()*vf()
                  - (
                        coefft0
                       *alpha.oldTime()()*rho.oldTime()()
                       *vf.oldTime()()*mesh().V0()
                      - coefft00
                       *alpha.oldTime().oldTime()()*rho.oldTime().oldTime()()
                       *vf.oldTime().oldTime()()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft
                   *alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField() -
                    (
                        coefft0
                       *alpha.oldTime().boundaryField()
                       *rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()

                      - coefft00
                       *alpha.oldTime().oldTime().boundaryField()
                       *rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );
    }
    else
    {
        return tmp<VolField<Type>>
        (
            VolField<Type>::New
            (
                ddtName,
                rDeltaT*
                (
                    coefft*alpha*rho*vf
                  - coefft0*alpha.oldTime()*rho.oldTime()*vf.oldTime()
                  + coefft00*alpha.oldTime().oldTime()
                   *rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDdtScheme<Type>::fvmDdt
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

    const scalar rDeltaT = 1.0/deltaT_();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDdtScheme<Type>::fvmDdt
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

    const scalar rDeltaT = 1.0/deltaT_();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT*rho.value())*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()*mesh().V0()
          - coefft00*vf.oldTime().oldTime().primitiveField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*rho.value()*
        (
            coefft0*vf.oldTime().primitiveField()
          - coefft00*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDdtScheme<Type>::fvmDdt
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

    const scalar rDeltaT = 1.0/deltaT_();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    fvm.diag() = (coefft*rDeltaT)*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().V0()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0*rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()
          - coefft00*rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDdtScheme<Type>::fvmDdt
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

    const scalar rDeltaT = 1.0/deltaT_();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    fvm.diag() =
        (coefft*rDeltaT)*alpha.primitiveField()*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*
        (
            coefft0
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().V0()

          - coefft00
           *alpha.oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = rDeltaT*mesh().V()*
        (
            coefft0
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()

          - coefft00
           *alpha.oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtUfCorr
(
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), (mesh().Sf() & Uf.oldTime()))
       *rDeltaT
       *(
            mesh().Sf()
          & (
                (coefft0*Uf.oldTime() - coefft00*Uf.oldTime().oldTime())
              - fvc::interpolate
                (
                    coefft0*U.oldTime() - coefft00*U.oldTime().oldTime()
                )
            )
        )
    );
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtPhiCorr
(
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
       *rDeltaT
       *(
            (coefft0*phi.oldTime() - coefft00*phi.oldTime().oldTime())
          - fvc::dotInterpolate
            (
                mesh().Sf(),
                coefft0*U.oldTime() - coefft00*U.oldTime().oldTime()
            )
        )
    );
}


template<class Type>
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& rhoUf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if
    (
        U.dimensions() == dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        const VolField<Type> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        const VolField<Type> rhoU00
        (
            rho.oldTime().oldTime()*U.oldTime().oldTime()
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + rhoUf.name() + ')',
            this->fvcDdtPhiCoeff
            (
                rhoU0,
                mesh().Sf() & rhoUf.oldTime(),
                rho.oldTime()
            )
           *rDeltaT
           *(
                mesh().Sf()
              & (
                    (
                        coefft0*rhoUf.oldTime()
                      - coefft00*rhoUf.oldTime().oldTime()
                    )
                  - fvc::interpolate(coefft0*rhoU0 - coefft00*rhoU00)
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && rhoUf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + rhoUf.name() + ')',
            this->fvcDdtPhiCoeff
            (
                U.oldTime(),
                mesh().Sf() & rhoUf.oldTime(),
                rho.oldTime()
            )
           *rDeltaT
           *(
                mesh().Sf()
              & (
                    (
                        coefft0*rhoUf.oldTime()
                      - coefft00*rhoUf.oldTime().oldTime()
                    )
                  - fvc::interpolate
                    (
                        coefft0*U.oldTime()
                      - coefft00*U.oldTime().oldTime()
                    )
                )
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
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimFlux
    )
    {
        const VolField<Type> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        const VolField<Type> rhoU00
        (
            rho.oldTime().oldTime()*U.oldTime().oldTime()
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), rho.oldTime())
           *rDeltaT
           *(
                (coefft0*phi.oldTime() - coefft00*phi.oldTime().oldTime())
              - fvc::dotInterpolate
                (
                    mesh().Sf(),
                    coefft0*rhoU0 - coefft00*rhoU00
                )
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimFlux
    )
    {
        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), rho.oldTime())
           *rDeltaT
           *(
                (coefft0*phi.oldTime() - coefft00*phi.oldTime().oldTime())
              - fvc::dotInterpolate
                (
                    mesh().Sf(),
                    coefft0*U.oldTime() - coefft00*U.oldTime().oldTime()
                )
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
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (U.dimensions() == dimVelocity && Uf.dimensions() == dimVelocity)
    {
        const volScalarField alphaRho0(alpha.oldTime()*rho.oldTime());
        const volScalarField alphaRho00
        (
            alpha.oldTime().oldTime()*rho.oldTime().oldTime()
        );

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
                    (
                        coefft0*fvc::interpolate(alphaRho0)*Uf.oldTime()
                      - coefft00*fvc::interpolate(alphaRho00)
                       *Uf.oldTime().oldTime()
                    )
                  - fvc::interpolate
                    (
                        coefft0*alphaRho0*U.oldTime()
                      - coefft00*alphaRho00*U.oldTime().oldTime()
                    )
                )
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
tmp<typename backwardDdtScheme<Type>::fluxFieldType>
backwardDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const fluxFieldType& phi
)
{
    const dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(U);

    const scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    const scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    const scalar coefft0  = coefft + coefft00;

    if (U.dimensions() == dimVelocity && phi.dimensions() == dimFlux)
    {
        const volScalarField alphaRho0(alpha.oldTime()*rho.oldTime());
        const volScalarField alphaRho00
        (
            alpha.oldTime().oldTime()*rho.oldTime().oldTime()
        );

        return fluxFieldType::New
        (
            "ddtCorr("
          + alpha.name() + rho.name() + ',' + U.name() + ',' + phi.name()
          + ')',
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
           *rDeltaT
           *(
                (
                    coefft0*fvc::interpolate(alphaRho0)*phi.oldTime()
                  - coefft00*fvc::interpolate(alphaRho00)
                   *phi.oldTime().oldTime()
                )
              - fvc::dotInterpolate
                (
                    mesh().Sf(),
                    coefft0*alphaRho0*U.oldTime()
                  - coefft00*alphaRho00*U.oldTime().oldTime()
                )
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
tmp<surfaceScalarField> backwardDdtScheme<Type>::meshPhi
(
    const VolField<Type>& vf
)
{
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    const scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    const scalar coefftn_0 = 1 + coefft0_00;

    return surfaceScalarField::New
    (
        mesh().phi().name(),
        coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime()
    );
}


template<class Type>
tmp<scalarField> backwardDdtScheme<Type>::meshPhi
(
    const VolField<Type>& vf,
    const label patchi
)
{
    const scalar deltaT = deltaT_();
    const scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    const scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    const scalar coefftn_0 = 1 + coefft0_00;

    return
    (
        coefftn_0*mesh().phi().boundaryField()[patchi]
      - coefft0_00*mesh().phi().oldTime().boundaryField()[patchi]
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
