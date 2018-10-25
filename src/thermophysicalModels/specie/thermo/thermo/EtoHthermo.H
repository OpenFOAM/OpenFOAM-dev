inline scalar Cp
(
    const scalar p,
    const scalar T
) const
{
    return Cv(p, T) + EquationOfState::CpMCv(p, T);
}

inline scalar Hs
(
    const scalar p,
    const scalar T
) const
{
    return Es(p, T) + p/EquationOfState::rho(p, T);
}

inline scalar Ha
(
    const scalar p,
    const scalar T
) const
{
    return Ea(p, T) + p/EquationOfState::rho(p, T);
}
