tmp<surfaceScalarField> tddtPhi1;
tmp<surfaceScalarField> tddtPhi2;

if (faceMomentum)
{
    #include "pUf/DDtU.H"
}
