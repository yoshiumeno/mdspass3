
#include <iostream>
#include <math.h>
#include <iomanip>
#include "PowellMinimization.h"
#include "ReaxFF.h"
#include "Matrix.h"
#include "Chromosome.h"

static const int    MAX_ITER_OUTER   = 100;
static const int    MAX_ITER_INNER   = 801;
static const int    MAX_ITER_BRACKET = 100;
static const int    MAX_ITER_BRENT   = 100;
static const double VERY_SMALL       = 1.0e-12;
static const double EPS              = 0.001;
static const double GOLD_RATIO       = 0.3819660;

bool PowellMinimization::Run(ReaxPotential& pot, const Input& input)
{
  if( PFMPI::s_master ) {
    cout << endl;
    cout << "=== Start Powell Minimization ===" << endl;
    cout << endl;
  }

  Chromosome chrom;
  chrom.ExtractFrom(pot);

  Vector vef1;
  double cost = ReaxFF::CalcCost(chrom, vef1);

  if( cost < VERY_SMALL ) {
    return true;
  }

  int nV = chrom.size();
  int nF = vef1.size();

  Matrix gamma(nF, nV);
  Matrix d(nV, nV);

  vector<double> vMin(nV);
  vector<double> vMax(nV);
  const vector<int>& viVar = pot.GetVariableIndex();
  for(int i = 0; i < nV; ++i) { 
    vMin[i] = pot[viVar[i]].m_dMin;
    vMax[i] = pot[viVar[i]].m_dMax;
  }

  int iOuter = 0;
  for( ; iOuter < MAX_ITER_OUTER; ++iOuter ) {

    if( !InitGamma(chrom, vef1, vMin, vMax, gamma, d) ) {
      return false;
    }

    Matrix les(nV, nV);
    Vector p(nV);
    SetupLinearEquationSystem(gamma, vef1, les, p);

    double cost_outer = cost;

    if( PFMPI::s_master ) {
      cout << "Iteration   Cost" << endl;
    }

    for( int iInner = 0; iInner < MAX_ITER_INNER; ++iInner ) {

      Vector q(nV);
      int info = les.SolveLinearEquation(p, q);

      if( 0 < info ) {
        cerr << "Failed to solve linear equation system. "
             << "Restart." << endl;
        break;
      }

      Vector delta = d * q;

      for( int i = 0; i < nV; ++i ) {
        delta[i] = fmax(delta[i], vMin[i] - chrom[i]);
        delta[i] = fmin(delta[i], vMax[i] - chrom[i]);
      }

      double cost_inner = cost;

      LinePoint min1;
      LinePoint min2;
      cost = LineMinimization(chrom, delta, vef1, min1, min2);
      
      int iMax = 0;
      double dMax = 0.0;
      for( int i = 0; i < nV; ++i ) {
        double d = fabs(p[i] * q[i]);
        if( dMax < d ) {
          iMax = i;
          dMax = d;
        }
      }

      if( !UpdateGamma(min1, min2, iMax, gamma, delta) ) {
        cerr << "Matrix gamma become singular, restart."
             << endl;
        break;
      }

      for(int i = 0; i < nV; ++i) {
        d[i][iMax] = delta[i];
      }

      SetupLinearEquationSystem(gamma, vef1, les, p);

      if( PFMPI::s_master ) {
        cout << setw( 7) << right << iInner;
        cout << setw(15) << right << fixed << setprecision(5) 
             << cost << endl;
      }

      if( nV <= iInner && 
          cost_inner - cost < input.m_dPowellPrecision ) {
        break;
      }
    };

    if( cost_outer - cost < input.m_dPowellPrecision ) {
      break;
    }
  }

  if( MAX_ITER_OUTER <= iOuter ) {
    Utils::Error("Reached maximum iteration.");
    return false;
  }
  else {
    if( PFMPI::s_master ) {
      cout << "Converged successfully." << endl;
    }
    return true;
  }
}

bool PowellMinimization::InitGamma(
  const Chromosome& chrom, const Vector& vef, 
  const vector<double>& vMin, const vector<double>& vMax,
  Matrix& gamma, Matrix& d)
{
  if( PFMPI::s_master ) {
    cout << "Initializing gamma ... ";
  }

  int nF = vef.size();
  int nV = chrom.size();

  d.SetIdentity();

  Vector vefTmp(nF, 0.0);

  for( int i = 0; i < nV; ++i ) {
    Chromosome chromTmp(chrom);
    
    double e = EPS * ( vMax[i] - vMin[i] );
    if( vMax[i] < chromTmp[i] + e ) e = -e;
    chromTmp[i] += e;

    double fb = ReaxFF::CalcCost(chromTmp, vefTmp);

    Vector vdiff = (vefTmp - vef) / e;
    double dNorm = vdiff.Norm();
    if( dNorm < VERY_SMALL ) {
      Utils::Error("Failed to initialize gamma.");
      return false;
    }
    vdiff.Normalize();

    for( int j = 0; j < nF; ++j ) {
      gamma[j][i] = vdiff[j];
    }

    d[i][i] /= dNorm;
  }

  if( PFMPI::s_master ) {
    cout << "done" << endl;
  }

  return true;
}

void PowellMinimization::SetupLinearEquationSystem(
  const Matrix& gamma, const Vector& vef, Matrix& les, Vector& p)
{
  Matrix gammaT = gamma.Trans();
  p   = - gammaT * vef;
  les = gammaT * gamma;
}

double PowellMinimization::LineMinimization(
  Chromosome& chrom, const Vector& delta, Vector& vef1,
  LinePoint& min1, LinePoint& min2)
{
  const double BX = 0.1;

  Chromosome chromB(chrom);
  chromB += BX * delta;

  LinePoint l, c, r;
  l.f = chrom.m_dCost;
  l.v = vef1;
  
  r.x = BX;
  r.f = ReaxFF::CalcCost(chromB, r.v);

  Bracket(chrom, delta, l, c, r);
  
  Brent(l, c, r, chrom, delta, min1, min2);

  chrom += min1.x * delta;
  vef1 = min1.v;

  return min1.f;
}

bool PowellMinimization::UpdateGamma(
  const LinePoint& min1, const LinePoint& min2, int j, 
  Matrix& gamma, Vector& delta)
{
  int m = gamma.m_nRow;
  
  Vector u = (min1.v - min2.v) / (min1.x - min2.x);
  double mu = u * min1.v / min1.f;
  Vector v = u - mu * min1.v;

  double d = v.Norm();
  if( d < VERY_SMALL ) {
    return false;
  }
  v.Normalize();

  for(int i = 0; i < m; ++i) {
    gamma[i][j] = v[i];
  }

  delta /= d;

  return true;
}

bool PowellMinimization::Bracket(
  const Chromosome& chrom0, const Vector& delta,
  LinePoint& l, LinePoint& c, LinePoint& r)
{
  Chromosome chrom(chrom0);

  if( l.f < r.f ) {
    c = l;
    l.x = (c.x - r.x) / GOLD_RATIO + r.x;
    chrom = chrom0 + l.x * delta;
    l.f = ReaxFF::CalcCost(chrom, l.v);
  }
  else {
    c = r;
    r.x = (c.x - l.x) / GOLD_RATIO + l.x;
    chrom = chrom0 + r.x * delta;
    r.f = ReaxFF::CalcCost(chrom, r.v);
  }
  
  int iIter = 0;
  for( ; iIter < MAX_ITER_BRACKET; ++iIter ) {
    if( l.f > c.f && c.f < r.f ) {
      return true;
    }

    if( l.f < c.f ) {
      r = c;
      c = l;
      l.x = (c.x - r.x) / GOLD_RATIO + r.x;
      chrom = chrom0 + l.x * delta;
      l.f = ReaxFF::CalcCost(chrom, l.v);
      continue;
    }
    if( c.f > r.f ){
      l = c;
      c = r;
      r.x = (c.x - l.x) / GOLD_RATIO + l.x;
      chrom = chrom0 + r.x * delta;
      r.f = ReaxFF::CalcCost(chrom, r.v);
      continue;
    }
    if( l.f == c.f ) {
      c = l;
      l.x = (c.x - r.x) / GOLD_RATIO + r.x;
      chrom = chrom0 + l.x * delta;
      l.f = ReaxFF::CalcCost(chrom, l.v);
    }
    else {
      c = r;
      r.x = (c.x - l.x) / GOLD_RATIO + l.x;
      chrom = chrom0 + r.x * delta;
      r.f = ReaxFF::CalcCost(chrom, r.v);
    }
  }

  Utils::Error("Failed to make bracket.");
  return false;
}

double PowellMinimization::Brent(
  const LinePoint& l, const LinePoint& c, const LinePoint& r, 
  const Chromosome& chrom0, const Vector& delta,
  LinePoint& min1, LinePoint& min2)
{
  const double TOL   = 1.0e-1;
  const int    ZEPS  = 1.0e-10;

  Chromosome chrom(chrom0);

  double a = fmin(l.x, r.x);
  double b = fmax(l.x, r.x);

  LinePoint x = c;
  LinePoint v = c;
  LinePoint w = c;

  double e = 0.0; // last movement

  for( int iter = 0; iter < MAX_ITER_BRENT; ++iter ) {
    double xm = 0.5 * ( a + b );
    double tol1 = TOL * fabs(x.x) + ZEPS;
    double tol2 = 2.0 * tol1;
    if( fabs(x.x - xm) <= (tol2 - 0.5 * (b - a)) ) {
      min1 = x;
      min2 = w;
      return min1.x;
    }
    double d = 0.0;
    if( tol1 < fabs(e) ) {
      double r = (x.x - w.x) * (x.f - v.f);
      double q = (x.x - v.x) * (x.f - w.f);
      double p = (x.x - v.x) * q - (x.x - w.x) * r;
      q = 2.0 * (q - r);
      if( 0.0 < q ) 
        p = -p;
      else
        q = -q;
      
      if( fabs(p) < fabs(0.5 * q * e) && 
          p > q * (a - x.x) && p < q * (b - x.x) ) {
        e = d;
        d = p / q;
        double ux = x.x + d;
        if( ux - a < tol2 || b - ux < tol2) {
          d = copysign(tol1, xm - x.x);
        }
      }
      else {
        e = x.x < xm ? b - x.x : a - x.x;
        d = GOLD_RATIO * e;
      }
    }
    else {
      e = x.x < xm ? b - x.x : a - x.x;
      d = GOLD_RATIO * e;
    }
    LinePoint u;
    u.x = x.x;
    u.x += tol1 <= fabs(d) ? d : copysign(tol1, d);
    chrom = chrom0 + u.x * delta;
    u.f = ReaxFF::CalcCost(chrom, u.v);

    if( x.f < u.f ) {
      if( u.x < x.x )
        a = u.x;
      else
        b = u.x;
      if( u.f <= w.f || w.x == x.x ) {
        v = w;
        w = u;
      }
      else if( u.f <= v.f || v.x == x.x || v.x == w.x ) {
        v = u;
      }
    }
    else {
      if( x.x <= u.x )
        a = x.x;
      else
        b = x.x;
      v = w;
      w = x;
      x = u;
    }
  }
  
  Utils::Error("Too many iterations in brent");
  
  min1 = x;
  return min1.f;
}
