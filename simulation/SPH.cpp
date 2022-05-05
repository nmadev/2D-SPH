#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <unistd.h>
#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <iostream>
#include <vector>
#include <math.h>   
#include <chrono>
using namespace std;

bool animate = false;
bool grid = false;
bool color = false;
int bin_size = 8;
float size_scale = 2.0f;
int windowX = 1024/size_scale;
int windowY = 1024/size_scale;
float max_height = 0.0f;
float eps = 16.0f;
float dt = 0.0005f;
float mass = 2.5f;
float rest_density = 100.0f;
float gas_constant = 2000.0f;
float gravity = -10000.0f;
float viscosity = 800.0f;

float kernel_dist = eps * eps;

// float kernel_poly6 = 4.0f / (M_PI * pow(eps, 8.0f));
// float kernel_spiky = -10.0f  / (M_PI * pow(eps, 5.0f));
// float kernel_viscosity = 40.0f / (M_PI * pow(eps, 5.0f));

float kernel_poly6 = 4.0f / (M_PI * pow(eps, 8.0f));
float kernel_spiky = 10.0f / (M_PI * pow(eps, 5.0f));
float kernel_viscosity = 10.0f / (9.0f * M_PI * pow(eps, 5.0f));


struct particle
{
    float x;
    float y;
    float vx;
    float vy;
    float density;
    float pressure;
    float forceX;
    float forceY;
};

vector<particle *> particles;

void stepEnvironmentVars()
{
    for (auto &n : particles)
    {
        float density_n = 0.0f;
        int d = 0;
        for (auto &m : particles)
        {
            float x = n->x - m->x;
            float y = n->y - m->y;
            float dist_sq = (x * x) + (y * y);

            if (dist_sq < kernel_dist)
            {
                d += 1;
                density_n += mass * pow(kernel_dist - dist_sq, 3.0f) * kernel_poly6;
            }
        }
        if (d == 0)
        {
            cout << d << endl;
        }
        n->density = density_n;
        n->pressure = gas_constant * (density_n - rest_density);
    }
}

void stepForces()
{
    for (auto &n : particles)
    {
        float pressureX = 0.0f; float pressureY = 0.0f;
        float viscosityX = 0.0f; float viscosityY = 0.0f;
        for (auto &m : particles)
        {
            float x = n->x - m->x;
            float y = n->y - m->y;
            if (x == 0 && y == 0)
            {
                continue;
            }
            float dist_sq = x * x + y * y;
            float dist = sqrt(dist_sq);
            float dx = x / dist;
            float dy = y / dist;

            if (dist < eps)
            {
                // cout << m->density << endl;
                pressureX += mass * -1.0 * dx * kernel_spiky * pow(eps - dist, 3.0f) * (n->pressure + m->pressure) / (2.0f * m->density);
                pressureY += mass * -1.0 * dy * kernel_spiky * pow(eps - dist, 3.0f) * (n->pressure + m->pressure) / (2.0f * m->density);
                float simplification_factor = pow(eps - dist, 2.0);
                viscosityX += mass * viscosity * (m->vx - n->vx) * kernel_viscosity * (simplification_factor) / (m->density);
                viscosityY += mass * viscosity * (m->vy - n->vy) * kernel_viscosity * (simplification_factor) / (m->density);
                // viscosityX += mass * viscosity * (m->vx - n->vx) * kernel_viscosity * (eps - dist) / (m->density);
                // viscosityY += mass * viscosity * (m->vy - n->vy) * kernel_viscosity * (eps - dist) / (m->density);
            }
            // cout << "PX: " << pressureX << " PY: " << pressureY << endl;
        }
        n->vx += (pressureX + viscosityX) * dt / n->density;
        n->vy += (gravity + pressureY + viscosityY) * dt / n->density;;
        // cout << "VX: " << n->vx << " VY:" << n->vy << endl;
    }
}

void imposeBoundaries()
{
    float reflection = 0.75;
    max_height = 0.0f;
    for (auto &n : particles)
    {
        n->x += n->vx * dt;
        n->y += n->vy * dt;

        if (n->x < eps)
        {
            // n->x = (eps - n->x) + eps;
            n->x = eps;
            n->vx *= -1.0f * reflection;
        }             
        if (n->x > windowX - eps)
        {
            // n->x = (windowX - eps) - (n->x - (windowX - eps)) ;
            n->x = windowX - eps;
            n->vx *=  -1.0f * reflection;
        }
        if (n->y < eps)
        {
            // n->y = (eps - n->y) + eps;
            n->y = eps;
            n->vy *=  -1.0f * reflection;
        }               
        if (n->y > windowY - eps)
        {
            // n->x = (windowX - eps) - (n->x - (windowX - eps)) ;
            n->y = windowY - eps;
            n->vy *=  -1.0f * reflection;
        }
        if (n->y > max_height)
        {
            max_height = n->y;
        }
    }
}

void initialConfig()
{
    for (float x = 2.0f * eps; x < windowX * 1.0f - 2.0f * eps; x += 0.5f * eps)
    {
        for (float y = 2.0f * eps; y < windowY * 0.5f; y += 0.5f * eps)
        {
            float randX = ((float) rand()/RAND_MAX) - 0.5f;
            float randY = ((float) rand()/RAND_MAX) - 0.5f;
            // randX = 0.0f; randY = 0.0f;
            float randScale = 4.0f;

            particle *p = new particle;
            p->x = x + randX * randScale;
            p->y = y + randY * randScale;
            if (p->y > max_height)
            {
                max_height = p->y;
            }
            p->vx = 0.0f;
            p->vy = 0.0f;
            p->density = 0.0;
            p->pressure = 0.0;
            particles.push_back(p);
        }
    }
    cout << "NUMBER OF PARTICLES: " << particles.size() << endl;
}


void initialGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_POINT_SMOOTH);
	glMatrixMode(GL_PROJECTION);
    glPointSize(8.0f);
}

void displayCallback()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity ();
	glViewport(0, 0, windowX, windowY);
    glOrtho(0, windowX, 0, windowY, 0, 1);
	// gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	// glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    if (grid)
    {
        int xBins = windowX / bin_size;
        int yBins = windowY / bin_size;
        float bins[xBins][yBins] = {0.0f};
        float max_counts = 0.0f;
        for (auto &n : particles)
        {
            int xBin = (int) (n->x / bin_size); 
            int yBin = (int) (n->y / bin_size);
            bins[xBin][yBin] += 1.0f;
            if (bins[xBin][yBin] > max_counts)
            {
                max_counts = bins[xBin][yBin];
            }
        }
        glBegin(GL_QUADS);
        for (int i = 0; i < xBins; i++)
        {
            for (int j = 0; j < yBins; j++)
            {
                float color_shade = bins[i][j] / max_counts;
                glColor4f(color_shade, color_shade, color_shade, 1.0f);
                glVertex2f(i       * bin_size, j       * bin_size);
                glVertex2f((i + 1) * bin_size, j       * bin_size);
                glVertex2f((i + 1) * bin_size, (j + 1) * bin_size);
                glVertex2f(i       * bin_size, (j + 1) * bin_size);
            }
        }
        glEnd();

        // float h = 1.0f / _xRes;
	    // glBegin(GL_QUADS);
		// for (int i = 1; i < _xRes; i++)  
        // {
        //     float x = (i - 0.5f) * h;
        //     for (int j = 1; j < _yRes; j++) 
        //     {
        //         float y = (j - 0.5f) * h;

        //         float density = _density(i,j);

        //         glColor3f(density, density, density);
        //         glVertex2f(x, y);
        //         glVertex2f(x + h, y);
        //         glVertex2f(x + h, y + h);
        //         glVertex2f(x, y + h);
        //     }
        // }
	    // glEnd();
    }
    else
    {
        glBegin(GL_POINTS);
        for (auto &n : particles)
        {
            float color_shade = 0.0f;
            if (color)
            {
                color_shade = n->y / max_height;
            }
            glColor4f(color_shade, color_shade, 1.0f, 1.0f);
            glVertex2f(n->x, n->y);
        }
        glEnd();
    }

	glutSwapBuffers();
}

void idleCallback()
{
    // cout << "IDLE " << particles.size() << endl;
    glutPostRedisplay();
    if (animate)
    {
        stepEnvironmentVars();
        stepForces();
        imposeBoundaries();
    }
    // sleep(3);
}

void keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'a':
		case 'A':
			animate = !animate;
			break;
		case 'c':
		case 'C':
      		color = !color;
			break;

        case 'g':
        case 'G':
            grid = !grid;
            break;

		case 'q':
		case 'Q':
			exit(0);
			break;
	}
}

int main(int argc, char ** argv)
{

    // initialize the window with dimensions, colors, and window name
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(windowX, windowY);
	glutCreateWindow(argv[0]);

	glutDisplayFunc(displayCallback);
    glutIdleFunc(idleCallback);
	glutKeyboardFunc(keyboardCallback);
    

    initialConfig();
    initialGL();

	// glutMouseFunc(mouseCallback);
	// glutMotionFunc(motionCallback);
	// glutReshapeFunc(reshapeCallback);

	glutMainLoop();
    return 0;
}
