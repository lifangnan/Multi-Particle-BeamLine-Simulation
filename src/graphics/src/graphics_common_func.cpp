#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstring>
#include <unistd.h>
#include "graphics_common_func.h"

GLfloat yellow[4] = {1.0, 1.0, 0.0, 1.0};
GLfloat white[4] = {1.0, 1.0, 1.0, 1.0};
GLfloat grey[4] = {0.2, 0.2, 0.2, 1.0};
GLfloat black[4] = {0, 0, 0, 1};
GLfloat red[4] = {1, 0, 0, 1};
GLfloat green[4] = {0, 1, 0, 1};
GLfloat blue[4] = {0, 0.4, 1, 1};
GLfloat skyblue[4] = {0.52, 0.74, 0.84, 1};
GLfloat magenta[4]= { 1.0, 0.0, 1.0 ,1};
GLfloat cyan[4]= { 0.0, 1.0, 1.0 ,1};
GLfloat orange[4]= { 1.0, 0.5, 0.0 ,1};

char* FileRead(const char* filename)
{
  FILE* in = fopen(filename, "rb");
  if (in == NULL) return NULL;

  int res_size = BUFSIZ;
  char* res = (char*)malloc(res_size);
  int nb_read_total = 0;

  while (!feof(in) && !ferror(in)) {
    if (nb_read_total + BUFSIZ > res_size) {
      if (res_size > 10 * 1024 * 1024) break;
      res_size = res_size * 2;
      res = (char*)realloc(res, res_size);
    }
    char* p_res = res + nb_read_total;
    nb_read_total += fread(p_res, 1, BUFSIZ, in);
  }

  fclose(in);
  res = (char*)realloc(res, nb_read_total + 1);
  res[nb_read_total] = '\0';
  return res;
}

void PrintLog(GLuint object)
{
  GLint log_length = 0;
  if (glIsShader(object))
    glGetShaderiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else if (glIsProgram(object))
    glGetProgramiv(object, GL_INFO_LOG_LENGTH, &log_length);
  else {
    fprintf(stderr, "printlog: Not a shader or a program\n");
    return;
  }

  char* log = (char*)malloc(log_length);

  if (glIsShader(object))
    glGetShaderInfoLog(object, log_length, NULL, log);
  else if (glIsProgram(object))
    glGetProgramInfoLog(object, log_length, NULL, log);

  fprintf(stderr, "%s", log);
  free(log);
}

GLint CreateShader(const char* filename, GLenum type)
{
  const GLchar* source = FileRead(filename);
  if (source == NULL) {
    fprintf(stderr, "Error opening %s: ", filename); perror("");
    return 0;
  }
  GLuint res = glCreateShader(type);
  glShaderSource(res, 1, &source, NULL);
  free((void*)source);

  glCompileShader(res);
  GLint compile_ok = GL_FALSE;
  glGetShaderiv(res, GL_COMPILE_STATUS, &compile_ok);
  if (compile_ok == GL_FALSE) {
    fprintf(stderr, "%s:", filename);
    PrintLog(res);
    glDeleteShader(res);
    return 0;
  }

  return res;
}

void PrintBitmapString(void* font, const char* s)
{
   if (s && strlen(s)) {
      while (*s) {
         glutBitmapCharacter(font, *s);
         s++;
      }
   }
}

uint CreateVBO(uint r_sz)
{
  GLuint vbo;
  // 1: num of buffer objects. 2: ID of the object
  glGenBuffers(1, &vbo); 
  // First arg is target, hint for vbo to decide mem location:sys, 
  // video card,or AGP. Second arg is ID.
  glBindBuffer(GL_ARRAY_BUFFER, vbo); 
  // Copy data to the buffer object. 
  // First arg is target, performance hint. 
  // Second is the size of data in bytes.
  // Thrid is the source data location. 
  // Fourth is the  usage flag, performance hint.
  glBufferData(GL_ARRAY_BUFFER, r_sz, 0, GL_DYNAMIC_DRAW);
  // Binding buffer with 0 swithes off vbo operation
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  return vbo;
}

float Lerp(float r_a, float r_b, float r_t)
{
  return r_a + r_t * (r_b - r_a);
}

void ColorRamp(float r_t, float* r_r)
{
  const int ncolors = 7;
///*
   float c[ncolors][3] = {
        { 1.0, 0.0, 0.0, },
        { 1.0, 0.5, 0.0, },
        { 1.0, 1.0, 0.0, },
        { 0.0, 1.0, 0.0, },
        { 0.0, 1.0, 1.0, },
        { 0.0, 0.0, 1.0, },
        { 1.0, 0.0, 1.0, }
      };
//*/
/*  
  float c[ncolors][3] = {
      { 1.0, 0.0, 1.0, },
      { 0.0, 0.0, 1.0, },
      { 0.0, 1.0, 1.0, },
      { 0.0, 1.0, 0.0, },
      { 1.0, 1.0, 0.0, },
      { 1.0, 0.5, 0.0, },
      { 1.0, 0.0, 0.0, }
  };
*/
  r_t = r_t * (ncolors-1);
  int i = (int) r_t;
  float u = r_t - floor(r_t);
  r_r[0] = Lerp(c[i][0], c[i + 1][0], u);
  r_r[1] = Lerp(c[i][1], c[i + 1][1], u);
  r_r[2] = Lerp(c[i][2], c[i + 1][2], u);
}

void ColorRampWhiteBackGround(float r_t, float* r_r, double r_max)
{
  if (r_t == 0.0)
  {
    r_r[0] = 0.0; r_r[1] = 0.0; r_r[2] = 0.0;
    return;
  }
  const int ncolors = 6;
  float c[ncolors][3] = {
//      { 1.0, 0.0, 0.5, },
      { 0.0, 0.2, 1.0, },
      { 0.0, 1.0, 1.0, },
      { 0.0, 1.0, 0.0, },
      { 1.0, 1.0, 0.0, },
      { 1.0, 0.5, 0.0, },
      { 1.0, 0.0, 0.0, }
  };
  r_t = sqrt(r_t) * (ncolors - 1);
  int i = (int) r_t;
  float u = r_t - floor(r_t);
  r_r[0] = Lerp(c[i][0], c[i + 1][0], u);
  r_r[1] = Lerp(c[i][1], c[i + 1][1], u);
  r_r[2] = Lerp(c[i][2], c[i + 1][2], u);
}

std::string GetProjectTopDir()
{
  const int PATH_MAX = 512;
  char buff[PATH_MAX];
  std::string path;
  if(getcwd(buff, PATH_MAX) != 0)
    path = std::string(buff);
  size_t found = path.find("hpsim");
  size_t found1 = path.find_first_of("/", found);
  return path.substr(0, found1);
}

/*
void ColorRampWhiteBackGround(float r_t, float* r_r, double r_max)
{
  float a = (1 - r_t) * 4.0;
  int x = floor(a);
  float y = a - x;
  switch(x)
  {
    case 0: 
      r_r[0] = 1.0; r_r[1] = y; r_r[2] = 0.0; break;
    case 1:
      r_r[0] = 1.0 - y; r_r[1] = 1.0; r_r[2] = 0.0; break;
    case 2:
      r_r[0] = 0.0; r_r[1] = 1.0; r_r[2] = y; break;
    case 3:
      r_r[0] = 0.0; r_r[1] = 1.0 - y; r_r[2] = 1.0; break;
    case 4:
      r_r[0] = 0.0; r_r[1] = 0.0; r_r[2] = 1.0; break;
  }
}
*/
