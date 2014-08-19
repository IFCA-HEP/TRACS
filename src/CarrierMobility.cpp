

JacoboniMobility::JacoboniMobility( char carrier_type, double T = 300)
{
  _T = T;
  if (carrier_type == 'e')
  {
    _mu0 = 1440. * std::pow(_T/300., -2.26);
    _vsat = 1.054e7  * std::pow(_T/300., -0.602);
    _beta = 1.075 * std::pow(_T/300., 0.220);
  }
  else if (carrier_type == 'h')
  {
    _mu0 = 474. * std::pow(_T/300., -2.619);
    _vsat = 0.940e7  * std::pow(_T/300., -0.226);
    _beta = 0.924 * std::pow(_T/300., 0.550);
  }
}

double JacoboniMobility::obtain_mobility(double e_field_mod)
{
  return _mu0/std::pow(1.0+std::pow(_mu0*e_field_mod/_vsat,_beta), 1.0/_beta)
}
