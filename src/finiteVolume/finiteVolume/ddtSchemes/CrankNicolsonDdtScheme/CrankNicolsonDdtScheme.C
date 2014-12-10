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

#include "CrankNicolsonDdtScheme.H"
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
template<class GeoField>
CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh
)
:
    GeoField(io, mesh),
    startTimeIndex_(-2) // This field is for a restart and thus correct so set
                        // the start time-index to correspond to a previous run
{
    // Set the time-index to the beginning of the run to ensure the field
    // is updated during the first time-step
    this->timeIndex() = mesh.time().startTimeIndex();
}


template<class Type>
template<class GeoField>
CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::DDt0Field
(
    const IOobject& io,
    const fvMesh& mesh,
    const dimensioned<typename GeoField::value_type>& dimType
)
:
    GeoField(io, mesh, dimType),
    startTimeIndex_(mesh.time().timeIndex())
{}


template<class Type>
template<class GeoField>
label CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
startTimeIndex() const
{
    return startTimeIndex_;
}


template<class Type>
template<class GeoField>
GeoField& CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
operator()()
{
    return *this;
}


template<class Type>
template<class GeoField>
void CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>::
operator=(const GeoField& gf)
{
    GeoField::operator=(gf);
}


template<class Type>
template<class GeoField>
CrankNicolsonDdtScheme<Type>::DDt0Field<GeoField>&
CrankNicolsonDdtScheme<Type>::ddt0_
(
    const word& name,
    const dimensionSet& dims
)
{
    if (!mesh().objectRegistry::template foundObject<GeoField>(name))
    {
        const Time& runTime = mesh().time();
        word startTimeName = runTime.timeName(runTime.startTime().value());

        if
        (
            (
                runTime.timeIndex() == runTime.startTimeIndex()
             || runTime.timeIndex() == runTime.startTimeIndex() + 1
            )
         && IOobject
            (
                name,
                startTimeName,
                mesh()
            ).headerOk()
        )
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        startTimeName,
                        mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh()
                )
            );
        }
        else
        {
            regIOobject::store
            (
                new DDt0Field<GeoField>
                (
                    IOobject
                    (
                        name,
                        mesh().time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh(),
                    dimensioned<typename GeoField::value_type>
                    (
                        "0",
                        dims/dimTime,
                        pTraits<typename GeoField::value_type>::zero
                    )
                )
            );
        }
    }

    DDt0Field<GeoField>& ddt0 = static_cast<DDt0Field<GeoField>&>
    (
        const_cast<GeoField&>
        (
            mesh().objectRegistry::template lookupObject<GeoField>(name)
        )
    );

    return ddt0;
}


template<class Type>
template<class GeoField>
bool CrankNicolsonDdtScheme<Type>::evaluate
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return ddt0.timeIndex() != mesh().time().timeIndex();
}

template<class Type>
template<class GeoField>
scalar CrankNicolsonDdtScheme<Type>::coef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (mesh().time().timeIndex() - ddt0.startTimeIndex() > 0)
    {
        return 1.0 + ocCoeff_;
    }
    else
    {
        return 1.0;
    }
}


template<class Type>
template<class GeoField>
scalar CrankNicolsonDdtScheme<Type>::coef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    if (mesh().time().timeIndex() - ddt0.startTimeIndex() > 1)
    {
        return 1.0 + ocCoeff_;
    }
    else
    {
        return 1.0;
    }
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicolsonDdtScheme<Type>::rDtCoef_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef_(ddt0)/mesh().time().deltaT();
}


template<class Type>
template<class GeoField>
dimensionedScalar CrankNicolsonDdtScheme<Type>::rDtCoef0_
(
    const DDt0Field<GeoField>& ddt0
) const
{
    return coef0_(ddt0)/mesh().time().deltaT0();
}


template<class Type>
template<class GeoField>
tmp<GeoField> CrankNicolsonDdtScheme<Type>::offCentre_
(
    const GeoField& ddt0
) const
{
    if (ocCoeff_ < 1.0)
    {
        return ocCoeff_*ddt0;
    }
    else
    {
        return ddt0;
    }
}


template<class Type>
const FieldField<fvPatchField, Type>& ff
(
    const FieldField<fvPatchField, Type>& bf
)
{
    return bf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + dt.name() + ')',
            dt.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + dt.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            ddtIOobject,
            mesh(),
            dimensioned<Type>
            (
                "0",
                dt.dimensions()/dimTime,
                pTraits<Type>::zero
            )
        )
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            dimensionedScalar rDtCoef0 = rDtCoef0_(ddt0);

            ddt0.dimensionedInternalField() =
            (
                (rDtCoef0*dt)*(mesh().V0() - mesh().V00())
              - mesh().V00()*offCentre_(ddt0.dimensionedInternalField())
            )/mesh().V0();
        }

        tdtdt().dimensionedInternalField() =
        (
            (rDtCoef*dt)*(mesh().V() - mesh().V0())
          - mesh().V0()*offCentre_(ddt0.dimensionedInternalField())
        )/mesh().V();
    }

    return tdtdt;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(vf - vf.oldTime()) - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*rho.value()*
                    (
                        mesh().V()*vf.internalField()
                      - mesh().V0()*vf.oldTime().internalField()
                    ) - mesh().V0()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*rho*(vf - vf.oldTime()) - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        mesh().V()*rho.internalField()*vf.internalField()
                      - mesh().V0()*rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - mesh().V00()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()*vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef*(rho*vf - rho.oldTime()*vf.oldTime())
              - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
CrankNicolsonDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
            alpha.dimensions()*rho.dimensions()*vf.dimensions()
        );

    IOobject ddtIOobject
    (
        "ddt(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    dimensionedScalar rDtCoef = rDtCoef_(ddt0);

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()
                   *alpha.oldTime().internalField()
                   *rho.oldTime().internalField()
                   *vf.oldTime().internalField()

                  - mesh().V00()
                   *alpha.oldTime().oldTime().internalField()
                   *rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                ) - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()

                  - alpha.oldTime().oldTime().boundaryField()
                   *rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDtCoef.dimensions()
               *alpha.dimensions()*rho.dimensions()*vf.dimensions(),
                (
                    rDtCoef.value()*
                    (
                        mesh().V()
                       *alpha.internalField()
                       *rho.internalField()
                       *vf.internalField()

                      - mesh().V0()
                       *alpha.oldTime().internalField()
                       *rho.oldTime().internalField()
                       *vf.oldTime().internalField()
                    ) - mesh().V00()*offCentre_(ddt0.internalField())
                )/mesh().V(),
                rDtCoef.value()*
                (
                    alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField()

                  - alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                ) - offCentre_(ff(ddt0.boundaryField()))
            )
        );
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                alpha.oldTime()
               *rho.oldTime()
               *vf.oldTime()

              - alpha.oldTime().oldTime()
               *rho.oldTime().oldTime()
               *vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDtCoef
               *(
                   alpha*rho*vf
                 - alpha.oldTime()*rho.oldTime()*vf.oldTime()
                )
              - offCentre_(ddt0())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + vf.name() + ')',
            vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.value()*mesh().V();

    vf.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*rho.value()*
                (
                    mesh().V0()*vf.oldTime().internalField()
                  - mesh().V00()*vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*rho.value()*
                (
                    vf.oldTime().boundaryField()
                  - vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*rho*(vf.oldTime() - vf.oldTime().oldTime())
                 - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*rho.value()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + rho.name() + ',' + vf.name() + ')',
            rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*rho.internalField()*mesh().V();

    vf.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()*rho.oldTime().internalField()
                   *vf.oldTime().internalField()
                  - mesh().V00()*rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                  - rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef*rho.oldTime().internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                rho.oldTime()*vf.oldTime()
              - rho.oldTime().oldTime()*vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef*rho.oldTime().internalField()*vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
CrankNicolsonDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + alpha.name() + ',' + rho.name() + ',' + vf.name() + ')',
            alpha.dimensions()*rho.dimensions()*vf.dimensions()
        );

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    scalar rDtCoef = rDtCoef_(ddt0).value();
    fvm.diag() = rDtCoef*alpha.internalField()*rho.internalField()*mesh().V();

    vf.oldTime().oldTime();
    alpha.oldTime().oldTime();
    rho.oldTime().oldTime();

    if (mesh().moving())
    {
        if (evaluate(ddt0))
        {
            scalar rDtCoef0 = rDtCoef0_(ddt0).value();

            ddt0.internalField() =
            (
                rDtCoef0*
                (
                    mesh().V0()
                   *alpha.oldTime().internalField()
                   *rho.oldTime().internalField()
                   *vf.oldTime().internalField()

                  - mesh().V00()
                   *alpha.oldTime().oldTime().internalField()
                   *rho.oldTime().oldTime().internalField()
                   *vf.oldTime().oldTime().internalField()
                )
              - mesh().V00()*offCentre_(ddt0.internalField())
            )/mesh().V0();

            ddt0.boundaryField() =
            (
                rDtCoef0*
                (
                    alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()

                  - alpha.oldTime().oldTime().boundaryField()
                   *rho.oldTime().oldTime().boundaryField()
                   *vf.oldTime().oldTime().boundaryField()
                )
              - offCentre_(ff(ddt0.boundaryField()))
            );
        }

        fvm.source() =
        (
            rDtCoef
           *alpha.oldTime().internalField()
           *rho.oldTime().internalField()
           *vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V0();
    }
    else
    {
        if (evaluate(ddt0))
        {
            ddt0 = rDtCoef0_(ddt0)*
            (
                alpha.oldTime()
               *rho.oldTime()
               *vf.oldTime()

              - alpha.oldTime().oldTime()
               *rho.oldTime().oldTime()
               *vf.oldTime().oldTime()
            ) - offCentre_(ddt0());
        }

        fvm.source() =
        (
            rDtCoef
           *alpha.oldTime().internalField()
           *rho.oldTime().internalField()
           *vf.oldTime().internalField()
          + offCentre_(ddt0.internalField())
        )*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            U.dimensions()
        );

    dimensionedScalar rDeltaT = rDtCoef_(ddt0);

    IOobject ddtIOobject
    (
        "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    fluxFieldType phiCorr
    (
        mesh().Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime()))
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff
            (
                U.oldTime(),
                mesh().Sf() & Uf.oldTime(),
                phiCorr
            )
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
        ddt0_<GeometricField<Type, fvPatchField, volMesh> >
        (
            "ddt0(" + U.name() + ')',
            U.dimensions()
        );

    dimensionedScalar rDeltaT = rDtCoef_(ddt0);

    IOobject ddtIOobject
    (
        "ddtCorr(" + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    fluxFieldType phiCorr
    (
        phi.oldTime() - (mesh().Sf() & fvc::interpolate(U.oldTime()))
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            ddtIOobject,
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh> >
            (
                "ddt0(" + rho.name() + ',' + U.name() + ')',
                U.dimensions()
            );

        dimensionedScalar rDeltaT = rDtCoef_(ddt0);

        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            mesh().Sf() & (Uf.oldTime() - fvc::interpolate(rhoU0))
        );

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff
                (
                    rhoU0,
                    mesh().Sf() & Uf.oldTime(),
                    phiCorr
                )
               *rDeltaT*phiCorr
            )
        );

        return ddtCorr;
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh> >
            (
                "ddt0(" + U.name() + ')',
                U.dimensions()
            );

        dimensionedScalar rDeltaT = rDtCoef_(ddt0);

        fluxFieldType phiCorr
        (
             mesh().Sf() & (Uf.oldTime() - fvc::interpolate(U.oldTime()))
        );

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    mesh().Sf() & Uf.oldTime(),
                    phiCorr
                )
               *rDeltaT*phiCorr
            )
        );

        return ddtCorr;
    }
    else
    {
        FatalErrorIn
        (
            "CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename CrankNicolsonDdtScheme<Type>::fluxFieldType>
CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh> >
            (
                "ddt0(" + rho.name() + ',' + U.name() + ')',
                U.dimensions()
            );

        dimensionedScalar rDeltaT = rDtCoef_(ddt0);

        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            phi.oldTime() - (mesh().Sf() & fvc::interpolate(rhoU0))
        );

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(rhoU0, phi.oldTime(), phiCorr)
               *rDeltaT*phiCorr
            )
        );

        return ddtCorr;
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        DDt0Field<GeometricField<Type, fvPatchField, volMesh> >& ddt0 =
            ddt0_<GeometricField<Type, fvPatchField, volMesh> >
            (
                "ddt0(" + U.name() + ')',
                U.dimensions()
            );

        dimensionedScalar rDeltaT = rDtCoef_(ddt0);

        fluxFieldType phiCorr
        (
            phi.oldTime() - (mesh().Sf() & fvc::interpolate(U.oldTime()))
        );

        tmp<fluxFieldType> ddtCorr
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
               *rDeltaT*phiCorr
            )
        );

        return ddtCorr;
    }
    else
    {
        FatalErrorIn
        (
            "CrankNicolsonDdtScheme<Type>::fvcDdtPhiCorr"
        )   << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    DDt0Field<surfaceScalarField>& meshPhi0 = ddt0_<surfaceScalarField>
    (
        "meshPhiCN_0",
        dimVolume
    );

    if (evaluate(meshPhi0))
    {
        meshPhi0 =
            coef0_(meshPhi0)*mesh().phi().oldTime() - offCentre_(meshPhi0());
    }

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            coef_(meshPhi0)*mesh().phi() - offCentre_(meshPhi0())
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
