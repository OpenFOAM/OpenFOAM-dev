inline scalar Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cv(p, T) + EquationOfState::CpMCv(p, T);
}

inline scalar hs
(
    const scalar p,
    const scalar T
) const
{
    return es(p, T) + p/EquationOfState::rho(p, T);
}

inline scalar ha
(
    const scalar p,
    const scalar T
) const
{
    return ea(p, T) + p/EquationOfState::rho(p, T);
}
