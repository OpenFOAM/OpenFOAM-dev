/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "CoEulerDdtScheme.H"
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
tmp<volScalarField> CoEulerDdtScheme<Type>::CorDeltaT() const
{
    const surfaceScalarField cofrDeltaT(CofrDeltaT());

    tmp<volScalarField> tcorDeltaT
    (
        volScalarField::New
        (
            "CorDeltaT",
            mesh(),
            dimensionedScalar(cofrDeltaT.dimensions(), 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    volScalarField& corDeltaT = tcorDeltaT.ref();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, facei)
    {
        corDeltaT[owner[facei]] =
            max(corDeltaT[owner[facei]], cofrDeltaT[facei]);

        corDeltaT[neighbour[facei]] =
            max(corDeltaT[neighbour[facei]], cofrDeltaT[facei]);
    }

    const surfaceScalarField::Boundary& cofrDeltaTbf =
        cofrDeltaT.boundaryField();

    forAll(cofrDeltaTbf, patchi)
    {
        const fvsPatchScalarField& pcofrDeltaT = cofrDeltaTbf[patchi];
        const fvPatch& p = pcofrDeltaT.patch();
        const labelUList& faceCells = p.patch().faceCells();

        forAll(pcofrDeltaT, patchFacei)
        {
            corDeltaT[faceCells[patchFacei]] = max
            (
                corDeltaT[faceCells[patchFacei]],
                pcofrDeltaT[patchFacei]
            );
        }
    }

    corDeltaT.correctBoundaryConditions();

    return tcorDeltaT;
}


template<class Type>
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::CofrDeltaT() const
{
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    const surfaceScalarField& phi =
        static_cast<const objectRegistry&>(mesh())
        .lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        surfaceScalarField Co
        (
            mesh().surfaceInterpolation::deltaCoeffs()
           *(mag(phi)/mesh().magSf())
           *deltaT
        );

        return max(Co/maxCo_, scalar(1))/deltaT;
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            static_cast<const objectRegistry&>(mesh())
           .lookupObject<volScalarField>(rhoName_).oldTime();

        surfaceScalarField Co
        (
            mesh().surfaceInterpolation::deltaCoeffs()
           *(mag(phi)/(fvc::interpolate(rho)*mesh().magSf()))
           *deltaT
        );

        return max(Co/maxCo_, scalar(1))/deltaT;
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);

        return tmp<surfaceScalarField>(nullptr);
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
CoEulerDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const volScalarField rDeltaT(CorDeltaT());

    const word ddtName("ddt("+dt.name()+')');

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            GeometricField<Type, fvPatchField, volMesh>::New
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
            rDeltaT.primitiveField()*dt.value()
           *(1.0 - mesh().Vsc0()/mesh().Vsc());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
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
tmp<GeometricField<Type, fvPatchField, volMesh>>
CoEulerDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    const word ddtName("ddt("+vf.name()+')');

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT()*
                (
                    vf()
                  - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
CoEulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT()*rho.value()*
                (
                    vf()
                  - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
CoEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    const word ddtName("ddt("+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT()*
                (
                    rho()*vf()
                  - rho.oldTime()()
                   *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
CoEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(CorDeltaT());

    const word ddtName("ddt("+alpha.name()+','+rho.name()+','+vf.name()+')');

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT()*
                (
                    alpha()
                   *rho()
                   *vf()

                  - alpha.oldTime()()
                   *rho.oldTime()()
                   *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*
                (
                    alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField()

                  - alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            GeometricField<Type, fvPatchField, volMesh>::New
            (
                ddtName,
                rDeltaT
               *(
                   alpha*rho*vf
                 - alpha.oldTime()*rho.oldTime()*vf.oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
CoEulerDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
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

    const scalarField rDeltaT(CorDeltaT()().primitiveField());

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
CoEulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
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

    const scalarField rDeltaT(CorDeltaT()().primitiveField());

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
CoEulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
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

    const scalarField rDeltaT(CorDeltaT()().primitiveField());

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
CoEulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
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

    const scalarField rDeltaT(CorDeltaT()().primitiveField());

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
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT()));

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
    (
        phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return fluxFieldType::New
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr)*rDeltaT*phiCorr
    );
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT()));

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
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(CorDeltaT()));

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr(phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0));

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
            this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr, rho.oldTime())
            *rDeltaT*phiCorr
        );
    }
    else if
    (
        U.dimensions() == dimDensity*dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr
        (
            phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return fluxFieldType::New
        (
            "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
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
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename CoEulerDdtScheme<Type>::fluxFieldType>
CoEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
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
        GeometricField<Type, fvPatchField, volMesh> rhoU0
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
tmp<surfaceScalarField> CoEulerDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return surfaceScalarField::New
    (
        "meshPhi",
        mesh(),
        dimensionedScalar(dimVolume/dimTime, 0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
