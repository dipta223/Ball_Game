#include "Settings.h"
#include "Scene.h"

#include <GL/freeglut.h>

#include <array>
#include <iostream>
#include <numbers>
#include <random>
#include <string>


// Textures
constexpr auto texture_count = 16;
extern std::array<unsigned int, texture_count> textures;


namespace Scene {
	// Global Variables

	// Key States
	static std::array<bool, 256> keyStates = { false };

	// Camera Variables
	static GLfloat cameraX = 0.0f;
	static GLfloat cameraY = 0.0f;
	static GLfloat cameraZ = 50.0f;

	static GLfloat lookX = 0.0f;
	static GLfloat lookY = 0.0f;
	static GLfloat lookZ = 0.0f;

	static GLfloat upX = 0.0f;
	static GLfloat upY = 1.0f;
	static GLfloat upZ = 0.0f;

	// Rotation around the axes
	static GLfloat angleX = 0.0f, angleY = 0.0f, angleZ = 0.0f;
	static GLfloat axisX = 0.0f, axisY = 0.0f, axisZ = 0.0f;
	static bool rotateX = false, rotateY = false, rotateZ = false;

	// Roll, Pitch, Yaw
	static GLdouble anglePitch = 90.0, angleYaw = 90.0, angleRoll = 0.0;

	// Lights Controls
	static bool light0 = true, light1 = false, light2 = false, light3 = true;
	static bool ambient = true, diffuse = true, specular = true;

	// Placements testing
	GLfloat x = 0.0f, y = 0.0f, z = 0.0f;

	// Ball Control
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 gen(rd()); // seed the generator
	std::uniform_int_distribution<> distr(-85, 85); // define the range
	
	GLfloat ballX = 0.0f, ballY = 85.0f, ballZ = 0.0f;
	GLfloat ballSpeedX = 0.0f, ballSpeedY = -0.1f, ballSpeedZ = 0.0f;

	// Basket Controls
	GLfloat basketX = 0.0f, basketY = -85.0f, basketZ = 0.0f;
	GLfloat basketSpeedX = 0.1f, basketSpeedY = 0.0f, basketSpeedZ = 0.0f;
	
	// Score
	int score = 0;

	// Random Variables
	int wired = 0;
	int animat = 0;
	const int nt = 40;				//number of slices along x-direction
	const int ntheta = 100;

	const int L = 5;
	GLfloat ctrlpoints[L + 1][3] =
	{
		{ 0.0, 0.0, 0.0},
		{ 0.0, 0.54, 0.0},
		{ 0.0, 1.01, 0.0},
		{ 0.54, 1.01, 0.0},
		{ 0.55, 0.55, 0.0},
		{ 0.55, 0.0, 0.0}
	};

	GLdouble modelview[16]; //var to hold the modelview info

	static void GetNormal3f(
		GLfloat x1, GLfloat y1, GLfloat z1,
		GLfloat x2, GLfloat y2, GLfloat z2,
		GLfloat x3, GLfloat y3, GLfloat z3)
	{
		GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

		Ux = x2 - x1;
		Uy = y2 - y1;
		Uz = z2 - z1;

		Vx = x3 - x1;
		Vy = y3 - y1;
		Vz = z3 - z1;

		Nx = Uy * Vz - Uz * Vy;
		Ny = Uz * Vx - Ux * Vz;
		Nz = Ux * Vy - Uy * Vx;

		glNormal3f(Nx, Ny, Nz);
	}

	static void SetNormal3f(
		GLfloat x1, GLfloat y1, GLfloat z1,
		GLfloat x2, GLfloat y2, GLfloat z2,
		GLfloat x3, GLfloat y3, GLfloat z3)
	{
		GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

		Ux = x2 - x1;
		Uy = y2 - y1;
		Uz = z2 - z1;

		Vx = x3 - x1;
		Vy = y3 - y1;
		Vz = z3 - z1;

		Nx = Uy * Vz - Uz * Vy;
		Ny = Uz * Vx - Ux * Vz;
		Nz = Ux * Vy - Uy * Vx;

		glNormal3f(-Nx, -Ny, -Nz);
	}

	static GLfloat quad_vertices[8][3] =
	{
		{-1.0, -1.0, -1.0},
		{ 1.0, -1.0, -1.0},
		{-1.0, -1.0,  1.0},
		{ 1.0, -1.0,  1.0},

		{-1.0,  1.0, -1.0},
		{ 1.0,  1.0, -1.0},
		{-1.0,  1.0,  1.0},
		{ 1.0,  1.0,  1.0}
	};

	static GLuint quad_indices[6][4] =
	{
		{0,2,3,1},
		{0,2,6,4},
		{2,3,7,6},
		{1,3,7,5},
		{1,5,4,0},
		{6,7,5,4}
	};

	// Draw Cube
	static void Cube(float R = 255, float G = 255, float B = 255, float alpha = 1)
	{
		float r = R / 255.0f, g = G / 255.0f, b = B / 255.0f;

		GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_ambient[] = { r, g, b, 1.0 };
		GLfloat mat_diffuse[] = { r, g, b, 1.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat mat_shininess[] = { 60 };

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

		/*if (em)
		{
			glMaterialfv(GL_FRONT, GL_EMISSION, matEm);
		}
		else
		{
			glMaterialfv(GL_FRONT, GL_EMISSION, noMat);
		}*/



		glBegin(GL_QUADS);
		for (GLint i = 0; i < 6; i++)
		{
			GetNormal3f(
				quad_vertices[quad_indices[i][0]][0], quad_vertices[quad_indices[i][0]][1], quad_vertices[quad_indices[i][0]][2],
				quad_vertices[quad_indices[i][1]][0], quad_vertices[quad_indices[i][1]][1], quad_vertices[quad_indices[i][1]][2],
				quad_vertices[quad_indices[i][2]][0], quad_vertices[quad_indices[i][2]][1], quad_vertices[quad_indices[i][2]][2]
			);

			glVertex3fv(&quad_vertices[quad_indices[i][0]][0]); glTexCoord2f(0, 0);
			glVertex3fv(&quad_vertices[quad_indices[i][1]][0]); glTexCoord2f(0, 1);
			glVertex3fv(&quad_vertices[quad_indices[i][2]][0]); glTexCoord2f(1, 1);
			glVertex3fv(&quad_vertices[quad_indices[i][3]][0]); glTexCoord2f(1, 0);
		}
		glEnd();
	}

	static long long nCr(int n, int r)
	{
		if (r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
		long long ans = 1;
		int i;

		for (i = 1; i <= r; i++)
		{
			ans *= n - r + i;
			ans /= i;
		}

		return ans;
	}

	static void BezierCurve(double t, float xy[2])
	{
		double y = 0;
		double x = 0;
		t = t > 1.0 ? 1.0 : t;
		for (int i = 0; i <= L; i++)
		{
			int ncr = nCr(L, i);
			double oneMinusTpow = pow(1 - t, double(L - i));
			double tPow = pow(t, double(i));
			double coef = oneMinusTpow * tPow * ncr;
			x += coef * ctrlpoints[i][0];
			y += coef * ctrlpoints[i][1];

		}
		xy[0] = float(x);
		xy[1] = float(y);
	}

	static void MinarCircleBezier()
	{
		int i, j;
		float x, y, z, r;				//current coordinates
		float x1, y1, z1, r1;			//next coordinates
		float theta;

		const float startx = 0, endx = ctrlpoints[L][0];
		//number of angular slices
		const float dx = (endx - startx) / nt;	//x step size
		const float dtheta = 2 * 3.1416 / ntheta;		//angular step size

		float t = 0;
		float dt = 1.0 / nt;
		float xy[2];
		BezierCurve(t, xy);
		x = xy[0];
		r = xy[1];
		//rotate about z-axis
		float p1x = 0.0, p1y = 0.0, p1z = 0.0, p2x = 0.0, p2y = 0.0, p2z = 0.0;
		for (i = 0; i < nt; ++i)  			//step through x
		{
			theta = 0;
			t += dt;
			BezierCurve(t, xy);
			x1 = xy[0];
			r1 = xy[1];

			//draw the surface composed of quadrilaterals by sweeping theta
			glBegin(GL_QUAD_STRIP);
			for (j = 0; j <= ntheta; ++j)
			{
				theta += dtheta;
				double cosa = cos(theta);
				double sina = sin(theta);
				y = r * cosa;
				y1 = r1 * cosa;	//current and next y
				z = r * sina;
				z1 = r1 * sina;	//current and next z

				//edge from point at x to point at next x
				glVertex3f(x, y, z);

				if (j > 0)
				{
					SetNormal3f(p1x, p1y, p1z, p2x, p2y, p2z, x, y, z);
				}
				else
				{
					p1x = x;
					p1y = y;
					p1z = z;
					p2x = x1;
					p2y = y1;
					p2z = z1;
				}

				glVertex3f(x1, y1, z1);
				//forms quad with next pair of points with incremented theta value
			}
			glEnd();
			x = x1;
			r = r1;
		} //for i
	}

	static void Curve()
	{
		const double t = glutGet(GLUT_ELAPSED_TIME) / 5000.0;
		const double a = t * 90.0;

		// glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (wired)
		{
			glPolygonMode(GL_FRONT, GL_LINE);
			glPolygonMode(GL_BACK, GL_LINE);

		}
		else
		{
			glPolygonMode(GL_FRONT, GL_FILL);
			glPolygonMode(GL_BACK, GL_FILL);
		}

		glPushMatrix();

		if (animat)
			glRotated(a, 0, 0, 1);

		//glRotatef(anglex, 1.0, 0.0, 0.0);
		//glRotatef(angley, 0.0, 1.0, 0.0);         	//rotate about y-axis
		//glRotatef(anglez, 0.0, 0.0, 1.0);

		glRotatef(90, 0.0, 0.0, 1.0);
		glTranslated(-3.5, 0, 0);
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //get the modelview info

		// void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);

		//  matColor(0.9,0.5,0.1,20);   // front face color
		// matColor(0.0,0.5,0.8,20,1);  // back face color

		glPushMatrix();
		// glRotatef(90,0,1,0);
		// glScalef(0.5,1,0.5);
		MinarCircleBezier();
		glPopMatrix();

		/*if (shcpt)
		{
			matColor(0.0, 0.0, 0.9, 20);
			showControlPoints();
		}*/

		glPopMatrix();
	}

	static void Circle3D(GLdouble radius)
	{
		GLUquadric* qobj = gluNewQuadric();
		gluQuadricTexture(qobj, GL_TRUE);

		glRotatef(270, 1, 0, 0);
		gluSphere(qobj, radius, 20, 20);
		gluDeleteQuadric(qobj);
	}

	static void Cylinder3D(GLdouble height, GLdouble rad, GLdouble rad2)
	{
		GLUquadric* qobj = gluNewQuadric();
		gluQuadricTexture(qobj, GL_TRUE);
		glRotatef(90, 1, 0, 0);

		gluCylinder(qobj, rad, rad2, height, 20, 20);
		gluDeleteQuadric(qobj);
	}

	// Texture Enum to Texture ID
	enum Textures
	{
		earth,
		water,
		cloud_sky,
		green_grass,
		black_road,
		red_brick,
		white_wall,
		wood,
		floor,
		roof_tile,
		tree,
		leaf
	};


	// Hall
	static GLfloat cube[8][3] =
	{
		{-1.0, -1.0, 1.0},
		{1.0, -1.0, 1.0},
		{1.0, 1.0, 1.0},
		{-1.0, 1.0, 1.0},


		{-1.0, -1.0, -1.0},
		{1.0, -1.0, -1.0},
		{1.0, 1.0, -1.0},
		{-1.0, 1.0, -1.0},
	};
	static GLubyte quadIndices[6][4] =
	{
		{0,1,2,3},
		{7,6,5,4},
		{2,6,7,3},

		{0,4,5,1},
		{2,1,5,6},
		{7,4,0,3},
	};

	static void DrawCube(GLint R = 255, GLint G = 255, GLint B = 255, GLboolean emission = false)
	{
		GLfloat r = R / 255.0f;
		GLfloat g = G / 255.0f;
		GLfloat b = B / 255.0f;

		GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_ambient[] = { r, g, b, 1.0 };
		GLfloat mat_diffuse[] = { r, g, b, 1.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat mat_shininess[] = { 60 };

		GLfloat mat_em[] = { r, g, b, 1.0 };

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

		if (emission) glMaterialfv(GL_FRONT, GL_EMISSION, mat_em);
		else glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);

		for (GLint i = 0; i < 6; i++)
			for (GLint i = 0; i < 6; i++)
			{
				glBegin(GL_QUADS);

				glVertex3fv(&cube[quadIndices[i][0]][0]);
				glTexCoord2f(1, 1);
				glVertex3fv(&cube[quadIndices[i][1]][0]);
				glTexCoord2f(1, 0);
				glVertex3fv(&cube[quadIndices[i][2]][0]);
				glTexCoord2f(0, 0);
				glVertex3fv(&cube[quadIndices[i][3]][0]);
				glTexCoord2f(0, 1);
				glEnd();
			}
	}

	static void Cylinder(GLint c1, GLint c2, GLint c3, GLboolean emission = false)
	{

		GLfloat r = c1 / 255.0;
		GLfloat g = c2 / 255.0;
		GLfloat b = c3 / 255.0;

		GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_ambient[] = { r, g, b, 1.0 };
		GLfloat mat_diffuse[] = { r, g, b, 1.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat mat_shininess[] = { 60 };

		GLfloat mat_em[] = { r,g,b,1.0 };

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

		if (emission) glMaterialfv(GL_FRONT, GL_EMISSION, mat_em);
		else glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);

		const double PI = 3.14159;

		/* top triangle */
		double i, resolution = 0.1;
		double height = 1;
		double radius = 0.5;

		glPushMatrix();
		glTranslatef(0, -0.5, 0);
		//top Circle
		glBegin(GL_TRIANGLE_FAN);
		glTexCoord2f(0.5, 0.5);
		glVertex3f(0, height, 0);  /* center */
		for (i = 2 * PI; i >= 0; i -= resolution)
		{
			glTexCoord2f(0.5f * cos(i) + 0.5f, 0.5f * sin(i) + 0.5f);
			glVertex3f(radius * cos(i), height, radius * sin(i));
		}
		/* close the loop back to 0 degrees */
		glTexCoord2f(0.5, 0.5);
		glVertex3f(radius, height, 0);
		glEnd();

		//bottom Circle
		glBegin(GL_TRIANGLE_FAN);
		glTexCoord2f(0.5, 0.5);
		glVertex3f(0, 0, 0);  /* center */
		for (i = 0; i <= 2 * PI; i += resolution)
		{
			glTexCoord2f(0.5f * cos(i) + 0.5f, 0.5f * sin(i) + 0.5f);
			glVertex3f(radius * cos(i), 0, radius * sin(i));
		}
		glEnd();

		//cylinder side
		glBegin(GL_QUAD_STRIP);
		for (i = 0; i <= 2 * PI; i += resolution)
		{
			const float tc = (i / (float)(2 * PI));
			glTexCoord2f(tc, 0.0);
			glVertex3f(radius * cos(i), 0, radius * sin(i));
			glTexCoord2f(tc, 1.0);
			glVertex3f(radius * cos(i), height, radius * sin(i));
		}
		/* close the loop back to zero degrees */
		glTexCoord2f(0.0, 0.0);
		glVertex3f(radius, 0, 0);
		glTexCoord2f(0.0, 1.0);
		glVertex3f(radius, height, 0);
		glEnd();

		glPopMatrix();
	}

	static void Wall()
	{
		GLfloat r = 0, g = 255, b = 0;
		
		// Back Wall
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[5]);
		glPushMatrix();
		glTranslatef(0, 0, -35);
		glScalef(100, 100, 1);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);

		// Left Wall
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[5]);
		glPushMatrix();
		glTranslatef(-100, 0, 0);
		glScalef(5, 100, 35);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);
		
		// Right Wall
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[5]);
		glPushMatrix();
		glTranslatef(100, 0, 0);
		glScalef(5, 100, 35);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);

		// Floor
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[7]);
		glPushMatrix();
		glTranslatef(0, -100, 0);
		glScalef(100, 5, 35);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);
		
		// Ceiling
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[7]);
		glPushMatrix();
		glTranslatef(0, 100, 0);
		glScalef(100, 5, 35);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);
	}

	static void Basket()
	{
		GLfloat r = 0, g = 255, b = 0;

		// Basket
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[6]);
		glPushMatrix();
		glTranslatef(basketX, -95, 25);
		glScalef(15, 5, 5);
		DrawCube();
		glPopMatrix();
		glDisable(GL_TEXTURE_2D);
	}
	
	static void Ball()
	{
		GLfloat r = 0, g = 255, b = 0;

		GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat mat_ambient[] = { r, g, b, 1.0 };
		GLfloat mat_diffuse[] = { r, g, b, 1.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat mat_shininess[] = { 100 };

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

		// Ball
		/*glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, textures[7]);*/
		glPushMatrix();
		glTranslatef(ballX, ballY, 25);
		glScalef(5, 5, 5);
		glutSolidSphere(2, 20, 20);
		glPopMatrix();
		//glDisable(GL_TEXTURE_2D);
	}

	static void DisplayText(std::string& str, int x, int y, int z)
	{
		GLfloat mat_ambient[] = { 1, 1, 1, 1.0 };
		GLfloat mat_diffuse[] = { 1, 1, 1, 1.0 };
		GLfloat mat_specular[] = { 1,1,1, 1.0 };
		GLfloat mat_shininess[] = { 10 };

		glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(1);
		//glColor3b(1,0,0);
		glPushMatrix();
		glTranslatef(x, y, z);
		glScalef(0.1, 0.1, 1);

		for (int i = 0; i < str.size(); i++)
		{
			glutStrokeCharacter(GLUT_STROKE_ROMAN, str[i]);
		}
		glPopMatrix();
	}
	
	static void Scene()
	{
		glPushMatrix();
		Wall();
		glPopMatrix();
		
		glPushMatrix();
		Basket();
		glPopMatrix();
		
		glPushMatrix();
		Ball();
		glPopMatrix();
		
		glPushMatrix();
		std::string str = "Score: " + std::to_string(score);
		DisplayText(str, 0, 0, 0);
		glPopMatrix();
	}

	static void Light()
	{
		//Light
		glEnable(GL_LIGHTING);

		GLfloat noLight[] = { 0, 0, 0, 1 };
		GLfloat lightAmb[] = { 0.5, 0.5, 0.5, 1 };
		GLfloat lightDif[] = { 1, 1, 1, 1 };
		GLfloat lightSpec[] = { 1, 1, 1, 1 };
		GLfloat light1Pos[] = { 0, -50, 0, 1 };
		GLfloat light4Pos[] = { 0,50,0, 1 };

		// GLfloat light1Pos[] = {90, 90, 90, 1};
		//  GLfloat light4Pos[] = {90, 90, -90, 1};
		GLfloat light2Pos[] = { 65, 30, -35, 1 }; //spot light
		GLfloat light3Pos[] = { -35, 40, -50, 1 }; //spot light  GLfloat light2Pos[] = {15, 40, -45, 1}; //spot light
		// GLfloat light3Pos[] = {15, 40, 45, 1

		glEnable(GL_LIGHT0);  //1
		glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDif);
		glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);
		glLightfv(GL_LIGHT0, GL_POSITION, light1Pos);


		glEnable(GL_LIGHT1); //2 spot light
		glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmb);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDif);
		glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpec);
		glLightfv(GL_LIGHT1, GL_POSITION, light2Pos);

		glEnable(GL_LIGHT2); //3 spot light
		glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmb);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDif);
		glLightfv(GL_LIGHT2, GL_SPECULAR, lightSpec);
		glLightfv(GL_LIGHT2, GL_POSITION, light3Pos);

		glEnable(GL_LIGHT3); //4
		glLightfv(GL_LIGHT3, GL_AMBIENT, lightAmb);
		glLightfv(GL_LIGHT3, GL_DIFFUSE, lightDif);
		glLightfv(GL_LIGHT3, GL_SPECULAR, lightSpec);
		glLightfv(GL_LIGHT3, GL_POSITION, light4Pos);

		//GLfloat spotDirection[] = {0, -18, 20, 1};

		//2    {90, 36, -120, 1};
		//3    {125, 36, 180, 1};

		GLfloat spotDirection[] = { 0, -1, 0, 1 };   //2
		GLfloat spotCutOff[] = { 60 };

		glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, spotDirection);
		glLightfv(GL_LIGHT2, GL_SPOT_CUTOFF, spotCutOff);

		GLfloat spotDirection2[] = { 0, -1, 0, 1 }; //3
		GLfloat spotCutOff2[] = { 60 };

		glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, spotDirection2);
		glLightfv(GL_LIGHT1, GL_SPOT_CUTOFF, spotCutOff2);


		if (light0)
		{
			glEnable(GL_LIGHT0);
		}
		else
		{
			glDisable(GL_LIGHT0);
		}

		if (light1)
		{
			glEnable(GL_LIGHT1);
		}
		else
		{
			glDisable(GL_LIGHT1);
		}

		if (light2)
		{
			glEnable(GL_LIGHT2);
		}
		else
		{
			glDisable(GL_LIGHT2);
		}
		if (light3)
		{
			glEnable(GL_LIGHT3);
		}
		else
		{
			glDisable(GL_LIGHT3);
		}


		if (ambient)
		{
			glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
			glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmb);
			glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmb);
			glLightfv(GL_LIGHT3, GL_AMBIENT, lightAmb);

		}
		else
		{
			glLightfv(GL_LIGHT0, GL_AMBIENT, noLight);
			glLightfv(GL_LIGHT1, GL_AMBIENT, noLight);
			glLightfv(GL_LIGHT2, GL_AMBIENT, noLight);
			glLightfv(GL_LIGHT3, GL_AMBIENT, noLight);
		}

		if (diffuse)
		{
			glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDif);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDif);
			glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDif);
			glLightfv(GL_LIGHT3, GL_DIFFUSE, lightDif);
		}
		else
		{
			glLightfv(GL_LIGHT0, GL_DIFFUSE, noLight);
			glLightfv(GL_LIGHT1, GL_DIFFUSE, noLight);
			glLightfv(GL_LIGHT2, GL_DIFFUSE, noLight);
			glLightfv(GL_LIGHT3, GL_DIFFUSE, noLight);
		}

		if (specular)
		{
			glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);
			glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpec);
			glLightfv(GL_LIGHT2, GL_SPECULAR, lightSpec);
			glLightfv(GL_LIGHT3, GL_SPECULAR, lightSpec);
		}
		else
		{
			glLightfv(GL_LIGHT0, GL_SPECULAR, noLight);
			glLightfv(GL_LIGHT1, GL_SPECULAR, noLight);
			glLightfv(GL_LIGHT2, GL_SPECULAR, noLight);
			glLightfv(GL_LIGHT3, GL_SPECULAR, noLight);
		}

	}

	static void Draw()
	{
		// Draw the scene
		glPushMatrix();
		Scene();
		glPopMatrix();
	}

	static void Keyboard()
	{
		if (Scene::keyStates['A'])
		{
			Scene::x -= 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}
		if (Scene::keyStates['D'])
		{
			Scene::x += 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}
		if (Scene::keyStates['W'])
		{
			Scene::z -= 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}
		if (Scene::keyStates['S'])
		{
			Scene::z += 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}
		if (Scene::keyStates['Q'])
		{
			Scene::y -= 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}
		if (Scene::keyStates['E'])
		{
			Scene::y += 0.5f;
#ifdef _DEBUG
			std::cout << "Translation: " << Scene::x << ", " << Scene::y << ", " << Scene::z << std::endl;
#endif // _DEBUG
		}

		// Basket Movement
		if (keyStates['['] && basketX > -80)
		{
			basketX -= basketSpeedX;
		}
		if (Scene::keyStates[']'] && basketX < 80)
		{
			basketX += basketSpeedX;
		}

		// Roll, Pitch, Yaw
		GLfloat x1 = lookX - cameraX;
		GLfloat z1 = lookZ - cameraZ;
		GLfloat r = sqrt(x1 * x1 + z1 * z1);

		GLfloat theta = 0.0f;
		if (x1 == 0)
		{
			if (z1 > 0)
				theta = 90;
			else if (z1 < 0)
				theta = -90;
		}
		else
			theta = atan(z1 / x1) * 180 / std::numbers::pi;

		if ((z1 > 0 && theta < 0) || (z1 < 0 && theta>0))
			theta += 180;
		else if (z1 < 0 && theta < 0)
			theta += 360;

		GLfloat turn_angle_step = 10;

		if (Scene::keyStates['u'])
		{
			theta -= turn_angle_step;
			theta = theta * std::numbers::pi / 180;
			lookX = r * cos(theta) + cameraX;
			lookZ = r * sin(theta) + cameraZ;
		}
		if (Scene::keyStates['U'])
		{
			theta += turn_angle_step;
			theta = theta * std::numbers::pi / 180;
			lookX = r * cos(theta) + cameraX;
			lookZ = r * sin(theta) + cameraZ;
		}

		if (Scene::keyStates['r'])
			lookY += 1.0f;

		if (Scene::keyStates['R'])
			lookY -= 1.0f;

		if (Scene::keyStates['p'])
		{
			anglePitch -= 1.0f;
			lookX = r * (cos(anglePitch * std::numbers::pi / 180.0)) * (cos(angleYaw * std::numbers::pi / 180.0));
			lookY = r * (sin(anglePitch * std::numbers::pi / 180.0));
			lookZ = r * (cos(anglePitch * std::numbers::pi / 180.0)) * (sin(angleYaw * std::numbers::pi / 180.0));
		}
		if (Scene::keyStates['P'])
		{
			anglePitch += 1.0f;
			lookX = 50.0 * (cos(anglePitch * std::numbers::pi / 180.0)) * (cos(angleYaw * std::numbers::pi / 180.0));
			lookY = 50.0 * (sin(anglePitch * std::numbers::pi / 180.0));
			lookZ = 50.0 * (cos(anglePitch * std::numbers::pi / 180.0)) * (sin(angleYaw * std::numbers::pi / 180.0));
		}

		GLfloat r1 = 1.0f;
		GLdouble dx = r1 * cos(theta * std::numbers::pi / 180);
		GLdouble dz = r1 * sin(theta * std::numbers::pi / 180);

		GLdouble dx_norm = r1 * cos((theta - 90) * std::numbers::pi / 180);
		GLdouble dz_norm = r1 * sin((theta - 90) * std::numbers::pi / 180);

		/*if (Scene::keyStates['j'])
		{
			cameraX += dx_norm * 3;
			cameraZ += dz_norm * 3;

			lookX += dx_norm * 3;
			lookZ += dz_norm * 3;
		}
		if (Scene::keyStates['i'])
		{
			cameraX += dx * 3;
			cameraZ += dz * 3;
			lookX += dx * 3;
			lookZ += dz * 3;
		}
		if (Scene::keyStates['k'])
		{
			cameraX -= dx * 3;
			cameraZ -= dz * 3;

			lookX -= dx * 3;
			lookZ -= dz * 3;
		}
		if (Scene::keyStates['l'])
		{
			cameraX -= dx_norm * 3;
			cameraZ -= dz_norm * 3;

			lookX -= dx_norm * 3;
			lookZ -= dz_norm * 3;
		}*/
	}

	static void Special()
	{
		if (keyStates[GLUT_KEY_LEFT])
		{
			Scene::cameraX -= 0.5f;
			Scene::lookX -= 0.5f;
		}
		if (keyStates[GLUT_KEY_RIGHT])
		{
			Scene::cameraX += 0.5f;
			Scene::lookX += 0.5f;
		}
		if (keyStates[GLUT_KEY_UP])
		{
			Scene::cameraY += 0.5f;
			Scene::lookY += 0.5f;
		}
		if (keyStates[GLUT_KEY_DOWN])
		{
			Scene::cameraY -= 0.5f;
			Scene::lookY -= 0.5f;
		}

		// Forward and Backward Movement
		if (keyStates[GLUT_KEY_PAGE_UP])
		{
			Scene::cameraZ -= 0.5f;
			Scene::lookZ -= 0.5f;
		}
		if (keyStates[GLUT_KEY_PAGE_DOWN])
		{
			Scene::cameraZ += 0.5f;
			Scene::lookZ += 0.5f;
		}
	}

	static void ResetBall()
	{
		ballX = distr(gen);
		ballY = 85;
	}
	
	static void Score()
	{
		score += 1;
		std::cout << "Score: " << score << std::endl;
		
		ResetBall();
	}
}

void Control::display()
{
	// Count the FPS
	static int frame = 0;
	static int time = 0;
	static int timebase = 0;
	char s[30]{};
	frame++;
	time = glutGet(GLUT_ELAPSED_TIME);
	if (time - timebase > 1000) {
		std::string fps = "FPS: " + std::to_string(frame * 1000.0 / (time - timebase));
		std::string title = Settings::WINDOW_TITLE + std::string(" - ") + fps;
		glutSetWindowTitle(title.c_str());
		//sprintf(s, "FPS:%4.2f", frame * 1000.0 / (time - timebase));
		timebase = time;
		frame = 0;
	}

	// Set up background color
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	// Clear the color and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set up the projection matrix
	gluOrtho2D(0, Settings::WINDOW_WIDTH, 0, Settings::WINDOW_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-10, 10, -10, 10, 5, 200);

	// Set up the modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(
		Scene::cameraX, Scene::cameraY, Scene::cameraZ,
		Scene::lookX, Scene::lookY, Scene::lookZ,
		Scene::upX, Scene::upY, Scene::upZ
	);

	glViewport(0, 0, Settings::WINDOW_WIDTH, Settings::WINDOW_HEIGHT);

	// Setting the rotation of the scene
	glRotatef(Scene::angleX, Scene::axisX, Scene::axisY, Scene::axisZ);
	glRotatef(Scene::angleY, Scene::axisX, Scene::axisY, Scene::axisZ);
	glRotatef(Scene::angleZ, Scene::axisX, Scene::axisY, Scene::axisZ);

	// Lighting
	Scene::Light();

	// Draw the scene
	Scene::Draw();

	// Keyboard input
	Scene::Keyboard();
	Scene::Special();

	// Swap the buffers
	glutSwapBuffers();

	// Redraw the scene
	glutPostRedisplay();
}

void Control::reshape(int width, int height) {}

void Control::keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
		// Scene Rotation (Positive Direction)
	case 'x':
		Scene::rotateX = !Scene::rotateX;
		Scene::axisX = 1.0f;
		break;
	case 'y':
		Scene::rotateY = !Scene::rotateY;
		Scene::axisY = 1.0f;
		break;
	case 'z':
		Scene::rotateZ = !Scene::rotateZ;
		Scene::axisZ = 1.0f;
		break;
		// Scene Rotation (Negative Direction)
	case 'X':
		Scene::rotateX = !Scene::rotateX;
		Scene::axisX = -1.0f;
		break;
	case 'Y':
		Scene::rotateY = !Scene::rotateY;
		Scene::axisY = -1.0f;
		break;
	case 'Z':
		Scene::rotateZ = !Scene::rotateZ;
		Scene::axisZ = -1.0f;
		break;

		// Light Control
	case '0':
		Scene::light0 = !Scene::light0;
		break;
	case '1':
		Scene::light1 = !Scene::light1;
		break;
	case '2':
		Scene::light2 = !Scene::light2;
		break;
	case '3':
		Scene::light3 = !Scene::light3;
		break;

		// Placement testing
	case 'A':
		Scene::keyStates['A'] = true;
		break;
	case 'D':
		Scene::keyStates['D'] = true;
		break;
	case 'W':
		Scene::keyStates['W'] = true;
		break;
	case 'S':
		Scene::keyStates['S'] = true;
		break;
	case 'Q':
		Scene::keyStates['Q'] = true;
		break;
	case 'E':
		Scene::keyStates['E'] = true;
		break;

		// Roll, Pitch, Yaw
	case 'r':
		Scene::keyStates['r'] = true;
		break;
	case 'p':
		Scene::keyStates['p'] = true;
		break;
	case 'u':
		Scene::keyStates['u'] = true;
		break;

		// Roll, Pitch, Yaw (Negative Direction)
	case 'R':
		Scene::keyStates['R'] = true;
		break;
	case 'P':
		Scene::keyStates['P'] = true;
		break;
	case 'U':
		Scene::keyStates['U'] = true;
		break;

		// Front-Back, Left-Right Movement
	/*case 'j':
		Scene::keyStates['j'] = true;
		break;
	case 'i':
		Scene::keyStates['i'] = true;
		break;
	case 'k':
		Scene::keyStates['k'] = true;
		break;
	case 'l':
		Scene::keyStates['l'] = true;
		break;*/
	
	case '[':
		Scene::keyStates['['] = true;
		break;
	case ']':
		Scene::keyStates[']'] = true;
		break;

		// Exit the program
	case 27:
		exit(EXIT_SUCCESS);
	default:
		break;
	}
}

void Control::keyboardup(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'A':
		Scene::keyStates['A'] = false;
		break;
	case 'D':
		Scene::keyStates['D'] = false;
		break;
	case 'W':
		Scene::keyStates['W'] = false;
		break;
	case 'S':
		Scene::keyStates['S'] = false;
		break;
	case 'Q':
		Scene::keyStates['Q'] = false;
		break;
	case 'E':
		Scene::keyStates['E'] = false;
		break;

	case 'r':
		Scene::keyStates['r'] = false;
		break;
	case 'p':
		Scene::keyStates['p'] = false;
		break;
	case 'u':
		Scene::keyStates['u'] = false;
		break;

	case 'R':
		Scene::keyStates['R'] = false;
		break;
	case 'P':
		Scene::keyStates['P'] = false;
		break;
	case 'U':
		Scene::keyStates['U'] = false;
		break;

		/*case 'j':
			Scene::keyStates['j'] = false;
			break;
		case 'i':
			Scene::keyStates['i'] = false;
			break;
		case 'k':
			Scene::keyStates['k'] = false;
			break;
		case 'l':
			Scene::keyStates['l'] = false;
			break;*/
		
	case '[':
		Scene::keyStates['['] = false;
		break;
	case ']':
		Scene::keyStates[']'] = false;
		break;

	default:
		break;
	}
}

void Control::special(int key, int x, int y)
{
	//std::cout << "Special Key Down Function called for " << key << std::endl;

	// Panning the camera
	switch (key)
	{
	case GLUT_KEY_LEFT:
		Scene::keyStates[GLUT_KEY_LEFT] = true;
		break;
	case GLUT_KEY_RIGHT:
		Scene::keyStates[GLUT_KEY_RIGHT] = true;
		break;
	case GLUT_KEY_UP:
		Scene::keyStates[GLUT_KEY_UP] = true;
		break;
	case GLUT_KEY_DOWN:
		Scene::keyStates[GLUT_KEY_DOWN] = true;
		break;

		// Zooming the camera
	case GLUT_KEY_PAGE_UP:
		Scene::keyStates[GLUT_KEY_PAGE_UP] = true;
		break;
	case GLUT_KEY_PAGE_DOWN:
		Scene::keyStates[GLUT_KEY_PAGE_DOWN] = true;
		break;
	default:
		break;
	}
}

void Control::specialup(int key, int x, int y)
{
	//std::cout << "Special Key Up Function called for " << key << std::endl;

	switch (key)
	{
	case GLUT_KEY_LEFT:
		Scene::keyStates[GLUT_KEY_LEFT] = false;
		break;
	case GLUT_KEY_RIGHT:
		Scene::keyStates[GLUT_KEY_RIGHT] = false;
		break;
	case GLUT_KEY_UP:
		Scene::keyStates[GLUT_KEY_UP] = false;
		break;
	case GLUT_KEY_DOWN:
		Scene::keyStates[GLUT_KEY_DOWN] = false;
		break;
		// Forward and Backward Movement
	case GLUT_KEY_PAGE_UP:
		Scene::keyStates[GLUT_KEY_PAGE_UP] = false;
		break;
	case GLUT_KEY_PAGE_DOWN:
		Scene::keyStates[GLUT_KEY_PAGE_DOWN] = false;
		break;
	default:
		break;
	}
}

void Control::mouse(int button, int state, int x, int y) {}

void Control::motion(int x, int y) {}

void Control::idle()
{
	if (Scene::rotateX)
	{
		Scene::angleX += 1.0f;
		if (Scene::angleX >= 360.0f)
			Scene::angleX = 0.0f;
	}
	if (Scene::rotateY)
	{
		if (Scene::angleY >= 360.0f)
			Scene::angleY = 0.0f;
		Scene::angleY += 1.0f;
	}
	if (Scene::rotateZ)
	{
		if (Scene::angleZ >= 360.0f)
			Scene::angleZ = 0.0f;
		Scene::angleZ += 1.0f;
	}

	// Ball Movement
	Scene::ballY += Scene::ballSpeedY;
	if (Scene::ballY <= -80 && abs(Scene::ballX - Scene::basketX) <= 10)
	{
		Scene::Score();
	}
	
	if (Scene::ballY <= -85) { 
		Scene::ResetBall();
	}

	// Redraw the scene
	glutPostRedisplay();
}
