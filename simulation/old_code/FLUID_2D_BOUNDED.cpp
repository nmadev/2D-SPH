///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D_BOUNDED.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D_BOUNDED::FLUID_2D_BOUNDED(int xRes, int yRes, float dt, float eps) :
  FLUID_2D(xRes, yRes, dt),
  _vorticityEps(eps),
  _residual(xRes, yRes),
  _direction(xRes, yRes),
  _q(xRes, yRes)
{
}

void FLUID_2D_BOUNDED::stomp_zero()
{
  for (int i = 0; i < _xRes; i++)
  {
    _xVelocity(0, i) = 0.0;
    _xVelocity(_xRes - 1, i) = 0.0;
    _yVelocity(0, i) = 0.0;
    _yVelocity(_xRes - 1, i) = 0.0;
  }
}

void FLUID_2D_BOUNDED::setNeumannBoundary(FIELD_2D& field)
{
  int N = _xRes;

  for (int i = 0; i < N; i++)
  {
    field(0,     i) = 0.0;
    field(i,     0) = 0.0;
    field(0, N - 1) = 0.0;
    field(N - 1, 0) = 0.0;
  }

  for (int i = 0; i < N; i++)
  {
    field(0,      i) = field(2,     i);
    field(N - 1,  i) = field(N - 2, i);
    field(i,      0) = field(i,     2);
    field(i,  N - 1) = field(i, N - 2);
  }
  // field(0,         0) = (field(0, 1) + field(1, 0)) / 2.0;
  // field(0,     N - 1) = (field(0, N - 2) + field(1, N - 1)) / 2.0;
  // field(N - 1,     0) = (field(N - 1, 1) + field(N - 2, 0)) / 2.0;
  // field(N - 1, N - 1) = (field(N - 1, N - 2) + field(N - 2, N - 1)) / 2.0;
  field(0,         0) = 0;
  field(0,     N - 1) = 0;
  field(N - 1,     0) = 0;
  field(N - 1, N - 1) = 0;
}

///////////////////////////////////////////////////////////////////////
// solve linear system with Gauss-Seidel iteration, with Dirichlet 
// boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::gaussSeidel(FIELD_2D& current, FIELD_2D& old)
{
  // IMPLEMENT ME
  for (int k = 0; k < 10; k++) 
  {
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
			  current(x,y) = (old(x,y) + current(x-1,y) + current(x+1,y) + current(x,y-1) + current(x,y+1)) * 0.25;

    setNeumannBoundary(current);
    // setPeriodicBoundary(current);
	}
}

///////////////////////////////////////////////////////////////////////
// advect field 'old' into 'current' using velocity field
// 'xVelocity' and 'yVelocity' and Dirichlet boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
{
  // IMPLEMENT ME
  int N = _xRes - 2;
	float dt0 = _dt * N;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      // trace backwards through the velocity field
		  float tempX = x - dt0 * xVelocity(x,y); 
      float tempY = y - dt0 * yVelocity(x,y);

      int   xInt    = (int)tempX;
      float xFrac = tempX - xInt;
      if (xInt > N) 
      {
        xInt = N;
      }
      if (xInt < 1) {
        xInt = 1;
      }
      tempX = xInt + xFrac;

      // wrap around for periodic boundary conditions in y direction
      int   yInt    = (int)tempY;
      float yFrac = tempY - yInt;
      if (yInt > N) 
      {
        yInt = N;
      }
      if (yInt < 1) {
        yInt = 1;
      }
      tempY = yInt + yFrac;

      // // wrap around for periodic boundary conditions in x direction
      // int   xInt    = (int)tempX;
      // float xFrac = tempX - xInt;
      // if (xInt > N) 
      //   xInt = xInt % N;
      // if (xInt < 1) {
      //   xInt = xInt % N;
      //   xInt += N - 1;
      // }
      // tempX = xInt + xFrac;

      // // wrap around for periodic boundary conditions in y direction
      // int   yInt    = (int)tempY;
      // float yFrac = tempY - yInt;
      // if (yInt > N) 
      //   yInt = yInt % N;
      // if (yInt < 1) {
      //   yInt = yInt % N;
      //   yInt += N - 1;
      // }
      // tempY = yInt + yFrac;

      // retrieve the coordinates of the grid cells to interpolate
      int x0 = (int) tempX; 
      int y0 = (int) tempY; 
      int x1 = x0 + 1;
      int y1 = y0 + 1;

      // compute the interpolation weights
		  float s1 = tempX - x0; 
      float s0 = 1 - s1; 
      float t1 = tempY - y0; 
      float t0 = 1 - t1;

      // compute the final interpolation
		  current(x,y) = s0 * (t0 * old(x0,y0) + t1 * old(x0, y1)) +
					           s1 * (t0 * old(x1,y0) + t1 * old(x1, y1));
    }
}

///////////////////////////////////////////////////////////////////////
// perform projection using Neumann boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::project()
{
  // IMPLEMENT ME
  int N = _xRes - 2;
  FIELD_2D& pressure = _xVelocityOld;
  FIELD_2D& divergence = _yVelocityOld;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      divergence(x,y) = -0.5f * (_xVelocity(x + 1,y)  - _xVelocity(x - 1,y) + 
                                 _yVelocity(x, y + 1) - _yVelocity(x, y - 1)) / N;
      pressure(x,y) = 0;
    }

  setNeumannBoundary(divergence);
  setNeumannBoundary(pressure);
	// gaussSeidel(pressure, divergence);
  // cout << "HERE" << endl;
  conjugateGradient(pressure, divergence);
	// gaussSeidel(pressure, divergence);
  // cout << "THERE" << endl;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      _xVelocity(x, y) -= 0.5f * N * (pressure(x + 1, y) - pressure(x - 1, y));
      _yVelocity(x, y) -= 0.5f * N * (pressure(x, y + 1) - pressure(x, y - 1));
    }
  // cout << "huh" << endl;

  setNeumannBoundary(_xVelocity);
  setNeumannBoundary(_yVelocity);
  // cout << "?????" << endl;
}

///////////////////////////////////////////////////////////////////////
// step density field using Neumann boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepDensity()
{
  // IMPLEMENT ME
  addSource(_density, _densityOld);
  swapFields(_density, _densityOld);
  stomp_zero();
  advect(_density, _densityOld, _xVelocity, _yVelocity);
  setNeumannBoundary(_density);
}

///////////////////////////////////////////////////////////////////////
// step velocity field using Neumann boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepVelocity()
{
  // cout << 1 << endl;
  addVorticity();
  // cout << 2 << endl;
  
  // cout << 3 << endl;
	addSource(_xVelocity, _xVelocityOld); 
  addSource(_yVelocity, _yVelocityOld);
	project();
  stomp_zero();

  // cout << 4 << endl;
	swapFields(_xVelocityOld, _xVelocity); 
  swapFields(_yVelocityOld, _yVelocity);
  
  // cout << 5 << endl;
	advect(_xVelocity, _xVelocityOld, _xVelocityOld, _yVelocityOld); 
  // cout << 5.5 << endl;
  setNeumannBoundary(_xVelocity);


  // cout << 6 << endl;
  stomp_zero();
  advect(_yVelocity, _yVelocityOld, _xVelocityOld, _yVelocityOld);
  setNeumannBoundary(_yVelocity);

  // cout << 7 << endl;
	project();
  
  // cout << 8 << endl;
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::conjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  int x, y, index;

  // i = 0
  int i = 0;
  int iterations = 10;

  // r = b - Ax
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
    {
      // if the cell is a variable
      float Acenter = 4.0f;
      
      _residual(x, y) = divergence(x, y) - (Acenter * pressure(x, y) -
        pressure(x - 1, y) - 
        pressure(x + 1, y) -
        pressure(x, y - 1) - 
        pressure(x, y + 1));
    }

  // d = r
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
      _direction(x, y) = _residual(x, y);

  // deltaNew = transpose(r) * r
  float deltaNew = 0.0f;
  for (y = 1; y < _yRes - 1; y++)
    for (x = 1; x < _xRes - 1; x++)
      deltaNew += _residual(x, y) * _residual(x, y);

  // delta0 = deltaNew
  // float delta0 = deltaNew;

  // While deltaNew > (eps^2) * delta0
  float eps  = 1e-08;
  float maxR = 2.0f * eps;

  while ((i < iterations) && (deltaNew > eps))
  {
    // q = Ad
    for (y = 1; y < _yRes - 1; y++, index += 2)
      for (x = 1; x < _xRes - 1; x++, index++)
      {
        // if the cell is a variable
        float Acenter = 4.0f;
        
        _q(x, y) = Acenter * _direction(x, y) -  
          _direction(x - 1, y) -
          _direction(x + 1, y) -
          _direction(x, y - 1) - 
          _direction(x, y + 1);
      }

    // alpha = deltaNew / (transpose(d) * q)
    float alpha = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        alpha += _direction(x, y) * _q(x, y);
    if (fabs(alpha) > 0.0f)
      alpha = deltaNew / alpha;

    // x = x + alpha * d
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        pressure(x, y) += alpha * _direction(x, y);

    // r = r - alpha * q
    maxR = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
      {
        _residual(x, y) -= alpha * _q(x, y);
        maxR = (_residual(x, y) > maxR) ? _residual(x, y) : maxR;
      }

    // deltaOld = deltaNew
    float deltaOld = deltaNew;

    // deltaNew = transpose(r) * r
    deltaNew = 0.0f;
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        deltaNew += _residual(x, y) * _residual(x, y);

    // beta = deltaNew / deltaOld
    float beta = deltaNew / deltaOld;

    // d = r + beta * d
    for (y = 1; y < _yRes - 1; y++)
      for (x = 1; x < _xRes - 1; x++)
        _direction(x, y) = _residual(x, y) + beta * _direction(x, y);
      
    setNeumannBoundary(pressure);

    // i = i + 1
    i++;
  } 
  cout << i << " iterations converged to " << maxR << endl;
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::preconditionedConjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  // IMPLEMENT ME
}

void FLUID_2D_BOUNDED::addVorticity()
{
  int x,y;
	if(_vorticityEps<=0.) return;

	// calculate vorticity
  // cout << "HERE1" << endl;
  for (y = 1; y < _yRes - 1; y++)
  {
    for (x = 1; x < _xRes - 1; x++)
    {
      float dx  = 1.0f / _xRes;
      float dy  = 1.0f / _yRes;
      _vorticity(x, y) = (_yVelocity(x + 1, y) - _yVelocity(x - 1, y)) * dy / 2.0
                            - (_xVelocity(x, y + 1) - _xVelocity(x, y - 1)) * dx / 2.0;
    }
  }
  float eps = _vorticityEps;

  for (y = 1; y < _yRes - 1; y++)
  {
    for (x = 1; x < _xRes - 1; x++)
    {
      float dx  = 1.0f / _xRes;
      float dy  = 1.0f / _yRes;
      float N[2];
      N[0] = (fabs(_vorticity(x + 1, y)) - fabs(_vorticity(x - 1, y))) * dx / 2.0; 
      N[1] = (fabs(_vorticity(x, y + 1)) - fabs(_vorticity(x, y - 1))) * dy / 2.0;

      float magnitude = sqrtf(N[0] * N[0] + N[1] * N[1]);
      if (magnitude > 0.0f)
      {
        magnitude = 1.0f / magnitude;
        N[0] *= magnitude;
        N[1] *= magnitude;
        _xVelocityOld(x, y) += (N[1] * _vorticity(x, y)) * _xRes * eps;
        _yVelocityOld(x, y) -= (N[0] * _vorticity(x, y)) * _yRes * eps;
      }
    }
  }
}