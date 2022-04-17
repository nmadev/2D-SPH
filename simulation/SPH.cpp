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

int windowX = 512;
int windowY = 512;
float eps = 4.0f;
float dt = 0.0001f;
float m = 1.0f;
float rest_density = 100.0f;
float gas_constant = 1000.0f;
float gravity = -10.0f;
// float viscosity = 10.0f / M_PI / pow(eps)
float kernel_dist = eps * eps;


struct particle
{
    float x;
    float y;
    float vx;
    float vy;
    float density;
    float pressure;
};

vector<particle *> particles;

// void stepEnvironmentVars()
// {
//     for (auto &n : particles)
//     {
//         float density_n = 0.0f
//         for auto (&m : particles)
//         {
//             float x = n->x - m->x;
//             float y = n->x - m->y;
//             float dist = sqrt(x * x + y * y);

//             if (dist < kernel_dist)
//             {
//                 density_n += m * 4.0f / (M_PI * pow(kernel_dist - dist, 4.0f));
//             }
//         }
//         n->density = density_n;
//         n->pressure = gas_constant * (density_n - rest_density);
//     }
// }

void stepForces()
{
    for (auto &n : particles)
    {
        n->vy += gravity * dt;
        // for (auto &m : particles)
        // {

        // }
    }
}

void imposeBoundaries()
{
    for (auto &n : particles)
    {
        n->x += n->vx * dt;
        n->y += n->vy * dt;

        if (n->x < eps)
        {
            n->x *= -1.0f;
            n->vx *= -0.9f;
        }             
        if (n->x > windowX - eps)
        {
            n->x *= -1.0f;
            n->vx *= -0.9f;
        }
        if (n->y < eps)
        {
            n->y *= -1.0f;
            n->vy *= -0.9f;
        }  
    }
}

void initialConfig()
{
    for (float x = eps; x < windowX * 1.0f; x += 2.0f * eps)
    {
        for (float y = eps; y < windowY * 0.5f; y += 2.0f * eps)
        {
            particle *p = new particle;
            p->x = x;
            p->y = y;
            p->vx = 0.0;
            p->vy = 0.0;
            p->density = 0.0;
            p->pressure = 0.0;
            particles.push_back(p);
        }
    }
}


void initialGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_POINT_SMOOTH);
	glMatrixMode(GL_PROJECTION);
    glPointSize(2.0f);
}

void displayCallback()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity ();
	glViewport(0, 0, windowX, windowY);
    glOrtho(0, windowX, 0, windowY, 0, 1);
	// gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	// glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // cout << "DISPLAY " << particles.size() << endl;

    glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
    glBegin(GL_POINTS);
    for (auto &n : particles)
    {
        glVertex2f(n->x, n->y);
    }
    glEnd();

	glutSwapBuffers();
}

void idleCallback()
{
    // cout << "IDLE " << particles.size() << endl;
    glutPostRedisplay();
    stepForces();
    imposeBoundaries();
    // sleep(3);
}

void keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		// case 'c':
		// case 'C':
      	// 	fluid->clear();
		// 	break;

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