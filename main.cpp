#include <GL/glew.h>
#include <GL/freeglut.h>
#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cstdlib>

#define PI 3.141092653589
#define DEG2RAD(deg) (deg * PI / 180)
#define RAD2DEG(rad) (rad * 180 / PI)
#define MAXSPEED 10
#define MINSPEED 0.5
#define MAXVEL 2
#define MINI 0.05
#define MAXI 2


char* loadFile(char, GLint);
void printShaderInfoLog(GLint);
void initShaders(void);
void handleKeypress1(unsigned char key, int x, int y);
void handleKeypress2(int key, int x, int y);
void handleMouseclick(int button, int state, int x, int y);
void motionCallBack(int x, int y);
void checkWin();

using namespace std;
using namespace glm;

GLuint p;
int level = 0;
bool first = true;
int score = 10;

struct Acc{
	double gx;
	double gy;
	double gz;
};

class Terrain{
	public:
		float width;
		float length;
		int tris;
		int divi;
		unsigned int VAO[205];
		unsigned int VBO[205];
	struct Acc findGravity(float x, float z)
	{
		float y;
		struct Acc gravity;
		double der = 0;
		gravity.gz = 0;
		gravity.gx = 0;
		gravity.gy = 0;
		if(level == 0)
		{
			y = 0.4*(x*x*x-3*x+z*z*z-3*z);
			//printf("x = %f y = %f z = %f\n",x,y,z );
			if(y > 0)
			{
			//	printf("sdfs\n");
				der = 1.0;
				//printf(" der = %lf, x = %f y = %f z = %f\n",der,x,y,z );
			}
			else if(y < 0)
			{
				//printf("sdfss\n");
				der = -1.0;
			}
			else
			{
				//printf("sdfsss\n");
				der = 0.0;
			}
			y = abs(y);
			//printf(" der = %lf, x = %f y = %f z = %f\n",der,x,y,z );
			gravity.gx = -0.4*(x*x*3 - 3)*der;
			gravity.gy = -0.4*der;
			gravity.gz = -0.4*(z*z*3 - 3)*der;
			//printf("der = %'f gx = %lf gy = %lf gz = %lf\n",der,gravity.gx,gravity.gy,gravity.gz);
			return gravity;

		}
		else if(level == 1)
			y = 0.4*fabs(cos(x*z)*(x*x-z*z));
		else if(level == 2)
			y = 0.4*fabs(sin(x*z)*(x*x-z*z) + x*x - z*x);
		else if(level == 3)
			y = 2*fabs(sin(x*z*z - x));
		else if(level == 4)
			y = 0.4*fabs(sin(x*z)/0.5+x*x + sin(x));


	}
	void setWidthHeight(float iw, float il)
	{
		width = iw;
		length = il;
		tris = 0;
	}
	float act_equation(float x, float z){
		if(level == 0)
			return 0.4*fabs((x*x*x-3*x+z*z*z-3*z));
		else if(level == 1)
			return 0.4*fabs(cos(x*z)*(x*x-z*z));
		else if(level == 2)
			return 0.4*fabs(sin(x*z)*(x*x-z*z) + x*x - z*x);
		else if(level == 3)
			return 2*fabs(sin(x*z*z - x));
		else if(level == 4)
			return 0.4*fabs(sin(x*z)/0.5+x*x + sin(x));
	}
	float equation(float i,float j){
		float x = i*width/divi;
		float z = ((int)(j/8))*length/divi;
		x = x - 2;
		z = z - 2;
		float y =  act_equation(x,z);
		return y;
	}
	float max(){
		return 3;
	}
	void makeTerrain(){
		GLfloat flatground[] = { 
									width,0,0,1,
									width,0,length,1,
									0,0,length,1,
									0,0,0,1
								};
		tris = 5;
		glGenVertexArrays(200,&VAO[0]);
		glBindVertexArray(VAO[0]);
		glGenBuffers(1,VBO);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, 16*sizeof(GLfloat), flatground, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		GLfloat trial[150][2*4*150] = {0};
		divi = 150;
		for (int i = 1; i <= divi; ++i)
		{
			for (int j = 0, k = 0; k < divi; j = j + 8, k++)
			{
				trial[i-1][j] = i*width/divi;
				trial[i-1][j + 1] = equation(i,j);
				trial[i-1][j + 2] = k*length/divi;
				trial[i-1][j + 3] = 1;
				trial[i-1][j + 4] = i*width/divi;
				trial[i-1][j + 5] = 0;
				trial[i-1][j + 6] = k*length/divi;
				trial[i-1][j + 7] = 1;
			}
			glBindVertexArray(VAO[i]);
			glGenBuffers(1,&VBO[i]);
			glBindBuffer(GL_ARRAY_BUFFER,VBO[i]);
			glBufferData(GL_ARRAY_BUFFER, 2*4*divi*sizeof(GLfloat), trial[i-1], GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);
		}
		glBindVertexArray(0);
	}
	void render(){
		float r = 0;
		float g = 0;
		float b = 0;
		glBindVertexArray(VAO[0]);
		glVertexAttrib3f(1,r/255,g/255,b/255);
		glDepthMask(GL_FALSE);
		glDrawArrays(GL_QUADS,0,4);
		glDepthMask(GL_TRUE);
		for (int i = 1; i <= divi; ++i)
		{
			//if(i > 60 && i <= 80 || i >= 40 && i <=45 || i > 50 && i < 55)
			//	continue;
			glBindVertexArray(VAO[i]);
			//237, 201, 175
			//150, 113, 23
			//244, 164, 96
					r = 237;
					g = 201;
					b = 175;
					float scale = equation(i,0);
					scale = scale/max();
					//scale = scale;
					r *= scale;
					g *= scale;
					b *= scale;
					float rk = 1,gk = 1,bk = 1;
					if(scale < 0.33)
						rk =2;
					else if(scale >= 0.33 && scale < 0.66)
						gk =2;
					else
						bk = 2;
			glVertexAttrib3f(1,rk*r/255,gk*g/255,bk*b/255);
			glDepthMask(GL_FALSE);
			glDrawArrays(GL_TRIANGLES,0,2*4*divi*divi);
			glDepthMask(GL_TRUE);
			glBindVertexArray(0);
		}	
	}
};

class Lattu{
	public:
		unsigned int VAO[10];
		unsigned int VBO[10];
		float x;
		float z;
		float y = 4;
		float r;
		float r2;
		float r3;
		float h1;
		float h2;
		float h3;
		float theta = 0.05;
		float phi;
		float velx = 0.5;
		float velz = 0.5;
		float ax = 0;
		float az = 0;
		float speed = 1;
		bool moving = false;
		float guideLineTheta = 0;
		float impulse;
		float x_axis_angle = 0;
		float y_axis_angle = 0;
		float z_axis_angle = 0;
	void set(float ix, float iz, float itheta, float iphi, Terrain terrain){
		//printf("received: %f %f\n",ix,iz);
		if(moving == true)
			checkWin();

		if(speed == 0)
		{
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0;
			theta = 0.5;
		}
		if(theta >= 0.8)
		{
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0;
			theta = 0.5;

		}
		if(ix > terrain.width)
		{
			ix = terrain.width;
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0;
			theta = 0;
		}
		if(iz > terrain.length)
		{
			iz = terrain.length;
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0.5;
		}
		if(ix < 0)
		{
			ix = 0;
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0;
			theta = 0.5;
		}
		if(iz < 0)
		{
			iz = 0;
			if(moving = true)
				score--;
			moving = false;
			velx = 0;
			velz = 0;
			theta = 0.5;
		}
		float yb = 0;
		x = ix;
		z = iz;
		yb = y;
		y = terrain.act_equation(x - terrain.width/2,z - terrain.length/2);
		yb = y - yb;
		if(fabs(velx) > MAXVEL )
			velx = MAXVEL;
		if(fabs(velz) > MAXVEL)
			velz = MAXVEL;

		if(fabs(velx) < MINI )
			velx = 0;
		if(fabs(velz) < MINI)
			velz = 0;

		if(moving == true)
		{
			// struct Acc gravity;
			// gravity = terrain.findGravity(x - terrain.width/2,z - terrain.length/2);
			// ax = (gravity.gx - x) * 0.02;
			// az = (gravity.gz - z) * 0.02;
			float olx,olz;
			olx = velx;
			olz = velz;
			if(yb > 0)
				ax = -0.08;
			else
				ax = 0.08;
			az = ax;
			//printf("%f %f %f\n",yb,ax,az );
			velx += ax;
			velz += az;
			if(olx*velx < 0)
				velx = 0;
			if(olz*velz < 0)
				velz = 0;
		}

		if(velx == 0 && velz == 0 && theta == 0.5 && moving == true)
		{
			score--;
			moving = false;
		}


		// iz = (int)iz/terrain.length;
		// ix = (int)ix/terrain.width;
		r = 0.13;
		r2 = 0.15;
		r3 = 0.07;
		h1 = 0.2;
		h2 = 0.3;
		h3 = 0.1;
		theta = itheta;
		phi = iphi;
		//axis(x,y,z);
		// if(moving == true)
		// 	printf("%f %f %f\n",x,y,z );
	}
	void transform(GLfloat *vertices, int size){
		// if(moving == false)
		// 	return;
		float tx, ty, tz;
		//printf("%d %f\n",moving, speed );
		for (int i = 0; i < size; i = i + 4)
		{
			tx = vertices[i] - x;
			ty = vertices[i + 1] - y;
			tz = vertices[i + 2] - z;
			vertices[i] = (tx*cos(phi) - tz*sin(phi))*cos(phi) + sin(phi)*(-ty*sin(theta) 
				+ cos(theta)*(tx*sin(phi) + tz*cos(phi))) + x;
			vertices[i + 1] = ty*cos(theta) + sin(theta)*(tx*sin(phi) + tz*cos(phi)) + y;
			vertices[i + 2] = -tx*cos(phi)*sin(phi) + tz*sin(phi)*sin(phi) - ty*sin(theta)*cos(phi) + 
			cos(theta)*cos(phi)*(tx*sin(phi) + tz*cos(phi)) + z;
			
		}
	}
	void display(float ix, float iz, float itheta, float iphi, Terrain terrain){
		//printf("sent: %f %f\n",ix,iz);
		set(ix,iz,itheta,iphi,terrain);
		GLfloat vertices[] = {
								x, y, z, 1,
								x + r*cos(DEG2RAD(360/6)), y + h1, z + r*sin(DEG2RAD(360/6)),1,
								x + r*cos(DEG2RAD(360/6 * 2)), y + h1, z + r*sin(DEG2RAD(360/6 * 2)),1,
								x + r*cos(DEG2RAD(360/6 * 3)), y + h1, z + r*sin(DEG2RAD(360/6 * 3)),1,
								x + r*cos(DEG2RAD(360/6 * 4)), y + h1, z + r*sin(DEG2RAD(360/6 * 4)),1,
								x + r*cos(DEG2RAD(360/6 * 5)), y + h1, z + r*sin(DEG2RAD(360/6 * 5)),1,
								x + r*cos(DEG2RAD(360/6 * 6)), y + h1, z + r*sin(DEG2RAD(360/6 * 6)),1,
								x, y, z, 1,
		};
		// Base of the hexagon thing
		GLfloat vertices2[] = {
								x, y + h1 , z, 1,
								x + r2*cos(DEG2RAD(360/6)), y + h1, z + r2*sin(DEG2RAD(360/6)),1,
								x + r2*cos(DEG2RAD(360/6 * 2)), y + h1, z + r2*sin(DEG2RAD(360/6 * 2)),1,
								x + r2*cos(DEG2RAD(360/6 * 3)), y + h1, z + r2*sin(DEG2RAD(360/6 * 3)),1,
								x + r2*cos(DEG2RAD(360/6 * 4)), y + h1, z + r2*sin(DEG2RAD(360/6 * 4)),1,
								x + r2*cos(DEG2RAD(360/6 * 5)), y + h1, z + r2*sin(DEG2RAD(360/6 * 5)),1,
								x + r2*cos(DEG2RAD(360/6 * 6)), y + h1, z + r2*sin(DEG2RAD(360/6 * 6)),1,
								x, y + h1 , z, 1,
		};
		// For the hexagon thing
		GLfloat vertices3[] = {
								x + r2*cos(DEG2RAD(360/6)), y + h1, z + r2*sin(DEG2RAD(360/6)),1,
								x + r2*cos(DEG2RAD(360/6)), y + h2, z + r2*sin(DEG2RAD(360/6)),1,
								x + r2*cos(DEG2RAD(360/6 * 2)), y + h2, z + r2*sin(DEG2RAD(360/6 * 2)),1,
								x + r2*cos(DEG2RAD(360/6 * 2)), y + h1, z + r2*sin(DEG2RAD(360/6 * 2)),1,

								x + r2*cos(DEG2RAD(360/6 * 2)), y + h2, z + r2*sin(DEG2RAD(360/6 * 2)),1,
								x + r2*cos(DEG2RAD(360/6 * 2)), y + h1, z + r2*sin(DEG2RAD(360/6 * 2)),1,
								x + r2*cos(DEG2RAD(360/6 * 3)), y + h1, z + r2*sin(DEG2RAD(360/6 * 3)),1,
								x + r2*cos(DEG2RAD(360/6 * 3)), y + h2, z + r2*sin(DEG2RAD(360/6 * 3)),1,

								x + r2*cos(DEG2RAD(360/6 * 3)), y + h1, z + r2*sin(DEG2RAD(360/6 * 3)),1,
								x + r2*cos(DEG2RAD(360/6 * 3)), y + h2, z + r2*sin(DEG2RAD(360/6 * 3)),1,
								x + r2*cos(DEG2RAD(360/6 * 4)), y + h2, z + r2*sin(DEG2RAD(360/6 * 4)),1,
								x + r2*cos(DEG2RAD(360/6 * 4)), y + h1, z + r2*sin(DEG2RAD(360/6 * 4)),1,

								x + r2*cos(DEG2RAD(360/6 * 4)), y + h2, z + r2*sin(DEG2RAD(360/6 * 4)),1,
								x + r2*cos(DEG2RAD(360/6 * 4)), y + h1, z + r2*sin(DEG2RAD(360/6 * 4)),1,
								x + r2*cos(DEG2RAD(360/6 * 5)), y + h1, z + r2*sin(DEG2RAD(360/6 * 5)),1,
								x + r2*cos(DEG2RAD(360/6 * 5)), y + h2, z + r2*sin(DEG2RAD(360/6 * 5)),1,

								x + r2*cos(DEG2RAD(360/6 * 5)), y + h1, z + r2*sin(DEG2RAD(360/6 * 5)),1,
								x + r2*cos(DEG2RAD(360/6 * 5)), y + h2, z + r2*sin(DEG2RAD(360/6 * 5)),1,
								x + r2*cos(DEG2RAD(360/6 * 6)), y + h2, z + r2*sin(DEG2RAD(360/6 * 6)),1,
								x + r2*cos(DEG2RAD(360/6 * 6)), y + h1, z + r2*sin(DEG2RAD(360/6 * 6)),1,

								x + r2*cos(DEG2RAD(360/6 * 6)), y + h2, z + r2*sin(DEG2RAD(360/6 * 6)),1,
								x + r2*cos(DEG2RAD(360/6 * 6)), y + h1, z + r2*sin(DEG2RAD(360/6 * 6)),1,
								x + r2*cos(DEG2RAD(360/6)), y + h1, z + r2*sin(DEG2RAD(360/6)),1,
								x + r2*cos(DEG2RAD(360/6)), y + h2, z + r2*sin(DEG2RAD(360/6)),1,
							
		};
		// For the handle
		GLfloat vertices4[] = {
							//	x, h3 + y + h1 , z, 1,
								x - r3*cos(DEG2RAD((360/6) * 1)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 1)),1,
								x - r3*cos(DEG2RAD((360/6) * 1)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 1)),1,
								x - r3*cos(DEG2RAD((360/6) * 2)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 2)),1,
								x - r3*cos(DEG2RAD((360/6) * 2)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 2)),1,

								x - r3*cos(DEG2RAD((360/6) * 2)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 2)),1,
								x - r3*cos(DEG2RAD((360/6) * 2)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 2)),1,
								x - r3*cos(DEG2RAD((360/6) * 3)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 3)),1,
								x - r3*cos(DEG2RAD((360/6) * 3)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 3)),1,

								x - r3*cos(DEG2RAD((360/6) * 3)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 3)),1,
								x - r3*cos(DEG2RAD((360/6) * 3)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 3)),1,
								x - r3*cos(DEG2RAD((360/6) * 4)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 4)),1,
								x - r3*cos(DEG2RAD((360/6) * 4)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 4)),1,

								x - r3*cos(DEG2RAD((360/6) * 4)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 4)),1,
								x - r3*cos(DEG2RAD((360/6) * 4)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 4)),1,
								x - r3*cos(DEG2RAD((360/6) * 5)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 5)),1,
								x - r3*cos(DEG2RAD((360/6) * 5)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 5)),1,

								x - r3*cos(DEG2RAD((360/6) * 5)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 5)),1,
								x - r3*cos(DEG2RAD((360/6) * 5)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 5)),1,
								x - r3*cos(DEG2RAD((360/6) * 6)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 6)),1,
								x - r3*cos(DEG2RAD((360/6) * 6)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 6)),1,

								x - r3*cos(DEG2RAD((360/6) * 6)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 6)),1,
								x - r3*cos(DEG2RAD((360/6) * 6)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 6)),1,
								x - r3*cos(DEG2RAD((360/6) * 1)), h3 + y + h1, z - r3*sin(DEG2RAD((360/6) * 1)),1,
								x - r3*cos(DEG2RAD((360/6) * 1)), h3 + y + h2, z - r3*sin(DEG2RAD((360/6) * 1)),1,
							//	x, h3 + y + h2 , z, 1,
		};

		transform(vertices, 32);
		transform(vertices2, 32);
		transform(vertices3, 6*4*4);
		transform(vertices4, 6*4*4);
		GLfloat colors1[] = {
								0,0,0,
								0,0.2,0,
								0,0.4,0,
								0,0.6,0,
								0,0.7,0,
								0,0.9,0,
								0,1.0,0,
								0,0,0,
		};
		GLfloat colors2[] = {
								0,0,0.1,
								0,0,0.2,
								0,0,0.3,
								0,0,0.4,
								0,0,0.4,
								0,0,0.3,
								0,0,0.2,
								0,0,0.1,
								0,0,0.1,
								0,0,0.2,
								0,0,0.3,
								0,0,0.4,
								0,0,0.4,
								0,0,0.3,
								0,0,0.2,
								0,0,0.1,
								0,0,0.1,
								0,0,0.2,
								0,0,0.3,
								0,0,0.4,
								0,0,0.4,
								0,0,0.3,
								0,0,0.2,
								0,0,0.1,						
		};
		GLfloat colors3[] = {
								0.1,0,0,
								0.2,0,0,
								0.3,0,0,
								0.4,0,0,
								0.4,0,0,
								0.3,0,0,
								0.2,0,0,
								0.1,0,0,
								0.1,0,0,
								0.2,0,0,
								0.3,0,0,
								0.4,0,0,
								0.4,0,0,
								0.3,0,0,
								0.2,0,0,
								0.1,0,0,
								0.1,0,0,
								0.2,0,0,
								0.3,0,0,
								0.4,0,0,
								0.4,0,0,
								0.3,0,0,
								0.2,0,0,
								0.1,0,0,
		};
		glGenVertexArrays(1,&VAO[0]);
		glBindVertexArray(VAO[0]);
		glGenBuffers(2,VBO);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, 8*4*sizeof(GLfloat), vertices, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, 8*3*sizeof(GLfloat), colors1, GL_STATIC_DRAW);
		glVertexAttribPointer((GLuint)1, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);

		glGenVertexArrays(1,&VAO[1]);
		glBindVertexArray(VAO[1]);
		glGenBuffers(1,&VBO[2]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[2]);
		glBufferData(GL_ARRAY_BUFFER, 8*4*sizeof(GLfloat), vertices2, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glGenVertexArrays(1,&VAO[2]);
		glBindVertexArray(VAO[2]);
		glGenBuffers(2,&VBO[3]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[3]);
		glBufferData(GL_ARRAY_BUFFER, 6*4*4*sizeof(GLfloat), vertices3, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, VBO[4]);
		glBufferData(GL_ARRAY_BUFFER, 6*4*3*sizeof(GLfloat), colors2, GL_STATIC_DRAW);
		glVertexAttribPointer((GLuint)1, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);

		glGenVertexArrays(1,&VAO[3]);
		glBindVertexArray(VAO[3]);
		glGenBuffers(2,&VBO[5]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[5]);
		glBufferData(GL_ARRAY_BUFFER, 6*4*4*sizeof(GLfloat), vertices4, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, VBO[6]);
		glBufferData(GL_ARRAY_BUFFER, 6*4*3*sizeof(GLfloat), colors3, GL_STATIC_DRAW);
		glVertexAttribPointer((GLuint)1, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(1);
	}
	void render()
	{
		glBindVertexArray(VAO[0]);
		glDepthMask(GL_FALSE);
		glDrawArrays(GL_TRIANGLE_FAN,0,8*4);
		glDepthMask(GL_TRUE);
		glBindVertexArray(0);

		glVertexAttrib3f(1,0,0,0.8);
		glBindVertexArray(VAO[1]);
		glDepthMask(GL_FALSE);
		glDrawArrays(GL_TRIANGLE_FAN,0,8*4);
		glDepthMask(GL_TRUE);
		glBindVertexArray(0);

		glBindVertexArray(VAO[2]);
		glDepthMask(GL_FALSE);
		glDrawArrays(GL_QUADS,0,6*4*4);
		glDepthMask(GL_TRUE);
		glBindVertexArray(0);

		glBindVertexArray(VAO[3]);
		glDepthMask(GL_FALSE);
		glDrawArrays(GL_QUADS,0,6*4*4);
		glDepthMask(GL_TRUE);
		glBindVertexArray(0);
	}
	void drawGuideLine()
	{
		//render();
		if(moving == false)
		{
			//float amp = sqrt(guideLineX*guideLineX + guideLineZ*guideLineZ);
			velx = impulse*cos(DEG2RAD(guideLineTheta));
			velz = impulse*sin(DEG2RAD(guideLineTheta));
			GLfloat guideLine[] = {
						x,y,z,1,
						x + velx, y, z + velz, 1,
					};
		
			//printf("Velx %f Velz %f\n",velx, velz );
			glGenVertexArrays(1,&VAO[4]);
			glBindVertexArray(VAO[4]);
			glGenBuffers(1,&VBO[7]);

			glBindBuffer(GL_ARRAY_BUFFER,VBO[7]);
			glBufferData(GL_ARRAY_BUFFER, 2*4*sizeof(GLfloat), guideLine, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,0,0.4,0);
			glBindVertexArray(VAO[4]);
			glDrawArrays(GL_LINES,0,2*4);
			glBindVertexArray(0);
			}
	}
	void move(Terrain terrain)
	{
		if(moving == true)
		{
			speed = speed/10;
			//printf("phi = %f theta = %f\n",phi, theta );
			display(x + velx ,z + velz ,theta + 0.0001*phi,phi + speed,terrain);
		}
	}
	void spin(Terrain terrain)
	{
		//printf("phi = %f theta = %f\n",phi, theta );
		display(x,z,theta + 0.0001*phi,phi+speed,terrain);
	}
};

class Target{
	public:
		unsigned int VAO[205];
		unsigned int VBO[205];
		float centerX;
		float centerY;
		float centerZ;
		float radius;
		float thickness;
		float theta;
	void set(Terrain terrain)
	{
		radius = 0.15;
		centerX = 4 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		centerZ = 4 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
		centerY = 0;
		thickness = 0.05;
		theta = 0;
		if(centerX > terrain.width)
			centerX = terrain.width;
		if(centerZ > terrain.length)
			centerZ = terrain.length;
		if(centerX < 0)
			centerX = 0;
		if(centerZ < 0)
			centerZ = 0;
		centerY = terrain.act_equation(centerX - terrain.width/2,centerZ - terrain.length/2) + radius;
		//printf("%f %f %f\n",centerX, centerY, centerZ );
	}
	void makeTarget(Terrain terrain)
	{
		set(terrain);
		float k = radius;
		GLfloat ring1[] = {
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,

								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
		};
		radius = (k*3)/4;
		GLfloat ring3[] = {
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,

								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
		};
		radius = (k*2)/4;
		GLfloat ring4[] = {
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,

								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
		};
		radius = (k*1)/4;
		GLfloat ring5[] = {
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								
								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,

								centerX + thickness, centerY , centerZ, 1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
		};
		// For the hecenterXagon thing
		radius = k;
		GLfloat ring2[] = {
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,

								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 2)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 2)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,

								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 3)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 3)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,

								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 4)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 4)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,

								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 5)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 5)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,

								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 6)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 6)),1,
								centerX, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
								centerX + thickness, centerY + radius*cos(DEG2RAD(theta + (360/6) * 1)), centerZ + radius*sin(DEG2RAD(theta + (360/6) * 1)),1,
							
		};
		glGenVertexArrays(1,&VAO[0]);
		glBindVertexArray(VAO[0]);
		glGenBuffers(1,VBO);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[0]);
		glBufferData(GL_ARRAY_BUFFER, 6*3*4*sizeof(GLfloat), ring1, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glGenVertexArrays(1,&VAO[1]);
		glBindVertexArray(VAO[1]);
		glGenBuffers(1,&VBO[1]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[1]);
		glBufferData(GL_ARRAY_BUFFER, 6*4*4*sizeof(GLfloat), ring2, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glGenVertexArrays(1,&VAO[2]);
		glBindVertexArray(VAO[2]);
		glGenBuffers(1,&VBO[2]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[2]);
		glBufferData(GL_ARRAY_BUFFER, 6*3*4*sizeof(GLfloat), ring3, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glGenVertexArrays(1,&VAO[3]);
		glBindVertexArray(VAO[3]);
		glGenBuffers(1,&VBO[3]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[3]);
		glBufferData(GL_ARRAY_BUFFER, 6*3*4*sizeof(GLfloat), ring4, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);

		glGenVertexArrays(1,&VAO[4]);
		glBindVertexArray(VAO[4]);
		glGenBuffers(1,&VBO[4]);

		glBindBuffer(GL_ARRAY_BUFFER,VBO[4]);
		glBufferData(GL_ARRAY_BUFFER, 6*3*4*sizeof(GLfloat), ring5, GL_STATIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
		glEnableVertexAttribArray(0);
	}
	void render()
	{
		glVertexAttrib3f(1,1,0,0);
		glBindVertexArray(VAO[0]);
		glDrawArrays(GL_TRIANGLES,0,6*3*4);
		glBindVertexArray(0);

		glVertexAttrib3f(1,1,1,1);
		glBindVertexArray(VAO[2]);
		glDrawArrays(GL_TRIANGLES,0,6*3*4);
		glBindVertexArray(0);

		glVertexAttrib3f(1,1,0,0);
		glBindVertexArray(VAO[3]);
		glDrawArrays(GL_TRIANGLES,0,6*3*4);
		glBindVertexArray(0);

		glVertexAttrib3f(1,1,1,1);
		glBindVertexArray(VAO[4]);
		glDrawArrays(GL_TRIANGLES,0,6*3*4);
		glBindVertexArray(0);

		glVertexAttrib3f(1,0.7,0,0);
		glBindVertexArray(VAO[1]);
		glDrawArrays(GL_QUADS,0,6*4*4);
		glBindVertexArray(0);

	}
};

class Camera{
	
	public:
		glm::mat4 Projection;
		glm::mat4 View;
		glm::mat4 Model;
		glm::mat4 MVP;	
		float camera_angle;
	public:
	void init()
	{
		camera_angle = 45;
		Projection = glm::perspective(camera_angle, 1.0f / 1.0f, 0.1f, 100.0f);
		View       = glm::lookAt(
		    glm::vec3(7,6,7), 
		    glm::vec3(0,0,0), 
		    glm::vec3(0,1,0) 
		);
		Model      = glm::mat4(1.0f);
		MVP        = Projection * View * Model;
		GLuint MatrixID = glGetUniformLocation(p, "MVP");
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	}
	void playerView()
	{
		Model      = glm::mat4(1.0f);
		View       = glm::lookAt(
		    glm::vec3(7,6,7), 
		    glm::vec3(0,0,0), 
		    glm::vec3(0,1,0) 
		);
		MVP        = Projection * View * Model;
		GLuint MatrixID = glGetUniformLocation(p, "MVP");
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	}
	void topView(Lattu top, Target target)
	{
		Model      = glm::mat4(1.0f);
		View       = glm::lookAt(
		    glm::vec3(top.x,top.y + top.h2 + top.h3,top.z), 
		    glm::vec3(target.centerX,target.centerY,target.centerZ), 
		    glm::vec3(0,1,0) 
		);
		MVP        = Projection * View * Model;
		GLuint MatrixID = glGetUniformLocation(p, "MVP");
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	}
	void overheadView()
	{
		Model      = glm::mat4(1.0f);
			View       = glm::lookAt(
			    glm::vec3(4,15,4), 
			    glm::vec3(0,0,0), 
			    glm::vec3(0,1,0) 
			);
			MVP        = Projection * View * Model;
			GLuint MatrixID = glGetUniformLocation(p, "MVP");
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);	
	}
	void helicopterCam(float zoom,float yaw, float pitch)
	{
		printf("zoom = %f\n",zoom);
		Model = glm::rotate(Model, yaw, glm::vec3(0, 1, 0));
		Model = glm::rotate(Model, pitch, glm::vec3(0, 0, 1));
		View       = glm::lookAt(
		    glm::vec3(2,10*zoom,2), 
		    glm::vec3(2,0,2), 
		    glm::vec3(-1,0,0) 
		);
		MVP        = Projection * View * Model;
		GLuint MatrixID = glGetUniformLocation(p, "MVP");
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	}
	void followCam(Lattu top)
	{
		Model      = glm::mat4(1.0f);
		View       = glm::lookAt(
		    glm::vec3(top.x - 1,top.y + top.h2 + top.h3,top.z - 1), 
		    glm::vec3(top.x,top.y,top.z), 
		    glm::vec3(0,1,0) 
		);
		MVP        = Projection * View * Model;
		GLuint MatrixID = glGetUniformLocation(p, "MVP");
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	}
};

Terrain terrain;
Lattu top;
Camera camera;
Target target;

class Game{
	public:
		int camera_count = 0;
		bool right_button_down = false;
		bool left_button_down = false;
		bool start = false;
		unsigned int scoreVAO[100];
		unsigned int scoreVBO[100];
		//int score = 10;
		float x = top.x;
		float y = top.y + top.h2 + top.h3 + 0.5;
		float z = top.z;
		float l = -0.1;
		float yl = 0.1;
		float offset = -0.15;
		bool game_over = false;
		float dragX;
		float dragZ;
		float zoomX = 1;
		bool update = false;
	public:
		void init()
		{
			x = -top.x;
			y = (top.y + top.h2 + top.h3 + 0.25);
			z = top.z;
		}
		void refresh()
		{
			level = static_cast <int> (rand());
			level = level%5;
			score = 10;
			terrain.makeTerrain();
			top.display(0.5,0.5,0,0,terrain);
			target.makeTarget(terrain);
			terrain.render();
		}
		void nextLevel()
		{
			level++;
			level = level%5;
			terrain.makeTerrain();
			top.display(0.5,0.5,0,0,terrain);
			target.makeTarget(terrain);
			terrain.render();

		}
		void updateCamera()
		{
			//printf("zoomXsd %f\n",zoomX);
			if(camera_count%5 == 0)
				camera.playerView();
			if(camera_count%5 == 1)
				camera.topView(top,target);
			if(camera_count%5 == 2)
				camera.overheadView();
			if(camera_count%5 == 3)
			{
				if(update == true)
				{
					printf("zoomBEFORE %f\n",zoomX);
					if(zoomX > 0)
						zoomX = (zoomX/1200);
					else
					{
						zoomX = fabs(zoomX/150);
						if(zoomX < 1)
							zoomX = 1.1;
					}
					if(zoomX == 0)
						zoomX = 1;
					printf("zoomAftere %f\n",zoomX);
					update = false;
				}
				camera.helicopterCam(zoomX,0,0);
			}
			if(camera_count%5 == 4)
				camera.followCam(top);
		}
		void handleKeypress1(unsigned char key, int x, int y)
		{
				if(key == 'c')
					camera_count++;
				if(key == 'r')
					refresh();
				if(key == 'm')
				{
					if(top.speed < MAXSPEED)
						top.speed += 0.1;	
				}
				if(key == 'l')
				{
					if(top.speed > 0)
						top.speed -= 0.1;	
				}
				if(key == 32)
					top.moving = true;
			    if(key == 27)
			    {
			    	printf("You quit!\n");
			    	exit(0);
			    }
		}
		void handleKeypress2(int key, int x, int y)
		{
			if (key == GLUT_KEY_RIGHT)
			{
				top.guideLineTheta += 5;
				
			}
			else if (key == GLUT_KEY_LEFT)
			{
				top.guideLineTheta -= 5;
				
			}
			else if (key == GLUT_KEY_UP)
			{
				if(top.impulse < MAXI)
					top.impulse += 0.1;	
			}
			else if (key == GLUT_KEY_DOWN)
			{
				if(top.impulse > 0.2)
					top.impulse -= 0.1;	
			}
		}
		void handleMouseclick(int button, int state, int x, int y)
		{
			if(camera_count%5 == 3)
			{
				switch(button){
				        case GLUT_LEFT_BUTTON:
				            if(state == GLUT_DOWN){
				            	left_button_down = true;
				            	dragX = x;
				            	dragZ = y;
				              
				            }
				            else if(state == GLUT_UP)
				            {
				            	left_button_down = false;
				            	dragX = x - dragX;
				            	dragZ = y - dragZ;
				            	//printf("dragX = %f dragZ = %f\n",0.4*dragX,0.4*dragZ );
				            	camera.helicopterCam(1,0.4*dragX,0.4*dragZ);
				            }
						case GLUT_RIGHT_BUTTON:
							if(state == GLUT_DOWN)
							{
								right_button_down = true;
								 zoomX = y;
								// dragZ = y;
							}
							else if(state = GLUT_UP)
							{
								zoomX = y - zoomX;
								if(right_button_down == true)
									update = true;
								right_button_down = false;
								// dragZ = y - dragZ;
								// printf("dragX = %f dragZ = %f\n",dragX,dragZ );
							}
				 		break;
				    }
			}
		}
		void motionCallBack(int x, int y)
		{
			if(right_button_down == false)
			{
				 glutPostRedisplay();
			}
		}
		void zero(float digit){
			GLfloat zero[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								x,y, offset*digit + z,1, // Lower left
								x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								/*x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,*/

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*6*4*sizeof(GLfloat), zero, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*6*4);
			glBindVertexArray(0);
		}
		void one(float digit){
			GLfloat one[] = {
								// x,y, offset*digit + z,1, // Bottom
								// x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Middle
								// x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Upper left
								// x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								// x,y + 2*yl, offset*digit + z,1, // Top
								// x,y + 2*yl, offset*digit + z + l,1,
			};
			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*2*4*sizeof(GLfloat), one, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*2*4);
			glBindVertexArray(0);
		}
		void two(float digit){
			GLfloat two[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								x,y, offset*digit + z,1, // Lower left
								x,y + yl, offset*digit + z,1,

								// x,y, offset*digit + z + l,1, // Lower right
								// x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Upper left
								// x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*5*4*sizeof(GLfloat), two, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*5*4);
			glBindVertexArray(0);
		}
		void three(float digit){
			GLfloat three[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Upper left
								// x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};


			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*5*4*sizeof(GLfloat), three, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*5*4);
			glBindVertexArray(0);
		}
		void four(float digit){
			GLfloat four[] = {
								// x,y, offset*digit + z,1, // Bottom
								// x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								// x,y + 2*yl, offset*digit + z,1, // Top
								// x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*4*4*sizeof(GLfloat), four, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*4*4);
			glBindVertexArray(0);
		}
		void five(float digit){
			GLfloat five[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								// x,y + yl, offset*digit + z + l,1, // Upper Right
								// x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*5*4*sizeof(GLfloat), five, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*5*4);
			glBindVertexArray(0);
		}
		void six(float digit){
			GLfloat six[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								x,y, offset*digit + z,1, // Lower left
								x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								// x,y + yl, offset*digit + z + l,1, // Upper Right
								// x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*6*4*sizeof(GLfloat), six, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*6*4);
			glBindVertexArray(0);
		}
		void seven(float digit){
			GLfloat seven[] = {
								// x,y, offset*digit + z,1, // Bottom
								// x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Middle
								// x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Upper left
								// x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*3*4*sizeof(GLfloat), seven, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*3*4);
			glBindVertexArray(0);
		}
		void eight(float digit){
			GLfloat eight[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								x,y, offset*digit + z,1, // Lower left
								x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*7*4*sizeof(GLfloat), eight, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*7*4);
			glBindVertexArray(0);
		}
		void nine(float digit){
			GLfloat nine[] = {
								// x,y, offset*digit + z,1, // Bottom
								// x,y, offset*digit + z + l,1,

								// x,y, offset*digit + z,1, // Lower left
								// x,y + yl, offset*digit + z,1,

								x,y, offset*digit + z + l,1, // Lower right
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Middle
								x,y + yl, offset*digit + z + l,1,

								x,y + yl, offset*digit + z,1, // Upper left
								x,y + 2*yl, offset*digit + z,1,

								x,y + yl, offset*digit + z + l,1, // Upper Right
								x,y + 2*yl, offset*digit + z + l,1,

								x,y + 2*yl, offset*digit + z,1, // Top
								x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*5*4*sizeof(GLfloat), nine, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*5*4);
			glBindVertexArray(0);
		}
		void starter(float digit){
			GLfloat nine[] = {
								x,y, offset*digit + z,1, // Bottom
								x,y, offset*digit + z + l,1,

								x,y, offset*digit + z,1, // Lower left
								x,y + yl, offset*digit + z,1,

								// x,y, offset*digit + z + l,1, // Lower right
								// x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Middle
								// x,y + yl, offset*digit + z + l,1,

								// x,y + yl, offset*digit + z,1, // Upper left
								// x,y + 2*yl, offset*digit + z,1,

								// x,y + yl, offset*digit + z + l,1, // Upper Right
								// x,y + 2*yl, offset*digit + z + l,1,

								// x,y + 2*yl, offset*digit + z,1, // Top
								// x,y + 2*yl, offset*digit + z + l,1,
			};

			glGenVertexArrays(1,&scoreVAO[0]);
			glBindVertexArray(scoreVAO[0]);
			glGenBuffers(1,scoreVBO);

			glBindBuffer(GL_ARRAY_BUFFER,scoreVBO[0]);
			glBufferData(GL_ARRAY_BUFFER, 2*2*4*sizeof(GLfloat), nine, GL_STATIC_DRAW);
			glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,0,0);
			glEnableVertexAttribArray(0);

			glVertexAttrib3f(1,1,0,0);
			glBindVertexArray(scoreVAO[0]);
			glDrawArrays(GL_LINES,0,2*2*4);
			glBindVertexArray(0);
		}
		void digitCaller(float digit, int number){
			if(number == 0)
				zero(digit);
			else if(number == 1)
				one(digit);
			else if(number == 2)
				two(digit);
			else if(number == 3)
				three(digit);
			else if(number == 4)
				four(digit);
			else if(number == 5)
				five(digit);
			else if(number == 6)
				six(digit);
			else if(number == 7)
				seven(digit);
			else if(number == 8)
				eight(digit);
			else if(number == 9)
				nine(digit);
		}
		void updateScore()
		{
			if(score >= 100)
				game_over = true;
			starter(0);
			int tempScore = score;
			digitCaller(2,tempScore%10);
			tempScore = tempScore/10;
			digitCaller(1,tempScore%10);	
		}
		float max(float x1, float x2)
		{
			if(x1 >= x2)
				return x1;
			return x2;
		}
		void checkWin()
		{
			//printf("top.x = %f top.y = %f top.x = %f\n",top.x,top.y,top.z );
			//printf("target.centerX = %f target.radius = %f target.centerZ = %f\n",top.x - target.centerX,target.radius,target.centerZ - top.z );
			if(fabs(top.x - target.centerX) <= (max(target.radius ,top.velx)) && 
				fabs(top.z - target.centerZ) <= (max(target.radius ,top.velz)))
			{
				score += 10;
				top.moving = false;
				nextLevel();

			}
		}

};

Game game;

void display(void)
{
	float r = 0;
	float g = 0;
	float b = 0;
	glClearColor(0.5f, 0.5f, 0.5f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	game.init();
	terrain.render();
	top.render();
	top.drawGuideLine();
	target.render();
	game.updateScore();
	game.updateCamera();
	glutSwapBuffers();
	top.speed = 5;
	// game.camera_count = 4;
	// top.spin(terrain);
	if(top.moving == true)
		top.move(terrain);
	if(game.game_over == true)
	{
		printf("GAMEOVER! Score: %d\n",score );
		exit(0);
	}
}

void idle()
{
	display();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glClearDepth (1.0f);
	glutInitWindowSize(600,600);
	glutCreateWindow("Terrain -- Lattu");
	glewInit();
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		cout << "glewInit failed, aborting." << endl;
		exit (1);
	}
	cout << "Status: Using GLEW " << glewGetString(GLEW_VERSION) << endl;
	cout << "OpenGL version " << glGetString(GL_VERSION) << " supported" << endl;
	initShaders();
	if(first)
	{
		printf("Welcome to Throwing Tops at Targets!\n");
		printf("To change the velocity at which you throw use the UP and DOWN arrow keys!\n");
		printf("To change the direction at which you throw use the LEFT and RIGHT arrow keys\n");
		printf("To change the speed/spin of the top use the keys 'm' and 'l' for MORE and LESS respectively\n");
		printf("To change the camera angle use the 'c' key.\t Note: The score may not be visible on all angles\n");
		printf("To make your move after setting the velocity,direction and spin, use the SPACE BAR\n");
		printf("To start a new game press 'r'\n" );
		printf("To exit the game press 'ESC'\n");
		terrain.setWidthHeight(4,4);
		terrain.makeTerrain();
		top.display(0.5,0.5,0,0,terrain);
		first = false;
		target.makeTarget(terrain);
	}
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(handleKeypress1);
    glutSpecialFunc(handleKeypress2);
    glutMouseFunc(handleMouseclick);
    glutMotionFunc(motionCallBack);
	glutPostRedisplay();	
	glutMainLoop();
	return 0;
}

char* loadFile(char *fname, GLint &fSize){
	ifstream::pos_type size;
	char * memblock;
	string text;

	// file read based on example in cplusplus.com tutorial
	ifstream file (fname, ios::in|ios::binary|ios::ate);
	if (file.is_open())
	{
		size = file.tellg();
		fSize = (GLuint) size;
		memblock = new char [size];
		file.seekg (0, ios::beg);
		file.read (memblock, size);
		file.close();
		cout << "file " << fname << " loaded" << endl;
		text.assign(memblock);
	}
	else
	{
		cout << "Unable to open file " << fname << endl;
		exit(1);
	}
	return memblock;
}

void printShaderInfoLog(GLint shader){
	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

	// should additionally check for OpenGL errors here

	if (infoLogLen > 0)
	{
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
		cout << "InfoLog:" << endl << infoLog << endl;
		delete [] infoLog;
	}

	// should additionally check for OpenGL errors here
}

void initShaders(void){
	GLuint f, v;

	char *vs,*fs;

	v = glCreateShader(GL_VERTEX_SHADER);
	f = glCreateShader(GL_FRAGMENT_SHADER);	

	// load shaders & get length of each
	GLint vlen;
	GLint flen;
	vs = loadFile("minimal.vert",vlen);
	fs = loadFile("minimal.frag",flen);
	
	const char * vv = vs;
	const char * ff = fs;

	glShaderSource(v, 1, &vv,&vlen);
	glShaderSource(f, 1, &ff,&flen);
	
	GLint compiled;

	glCompileShader(v);
	glGetShaderiv(v, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		cout << "Vertex shader not compiled." << endl;
		printShaderInfoLog(v);
	} 

	glCompileShader(f);
	glGetShaderiv(f, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		cout << "Fragment shader not compiled." << endl;
		printShaderInfoLog(f);
	} 
	
	p = glCreateProgram();

	glBindAttribLocation(p,0, "in_Position");
	glBindAttribLocation(p,1, "in_Color");
		
	glAttachShader(p,v);
	glAttachShader(p,f);
	
	glLinkProgram(p);
	glUseProgram(p);

	delete [] vs; // dont forget to free allocated memory
	delete [] fs; // we allocated this in the loadFile function...

	camera.init();
	//GLuint MatrixID = glGetUniformLocation(p, "MVP");
	//glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &camera.MVP[0][0]);
}

void handleKeypress1(unsigned char key, int x, int y){
	game.handleKeypress1(key,x,y);
}

void handleKeypress2(int key, int x, int y){
	game.handleKeypress2(key,x,y);
}

void handleMouseclick(int button, int state, int x, int y){
	game.handleMouseclick(button,state,x,y);
}

void motionCallBack(int x, int y){
	game.motionCallBack(x,y);
}

void checkWin()
{
	game.checkWin();
}

/*// DOES NOT WORK
	void axis(float x, float y, float z)
	{
		// equation is (x*x*x-3*x+z*z*z-3*z)
		// x normal is 3*x*x - 3
		// y normal is 0
		// z normal is 3*z*z - 3
		float xn = 0.4*fabs(3*x*x - 3);
		float yn = 0;
		float zn = 0.4*fabs(3*z*z - 3);
		xn = xn - x;
		yn = yn - y;
		zn = zn - z;
		float amp = sqrt(xn*xn + yn*yn + zn*zn);
		x_axis_angle = acos(xn/amp);
		y_axis_angle = acos(yn/amp);
		z_axis_angle = acos(zn/amp);
		// printf("%f %f %f\n",xn, yn, zn );
		// printf("%f %f %f\n",x_axis_angle, y_axis_angle, z_axis_angle );
	}*/