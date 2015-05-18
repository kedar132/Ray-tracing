//  Created by Hindupur Kedar on 11/18/14.
//  Copyright (c) 2014 Hindupur Kedar. All rights reserved.
//

#include <OpenGL/gl.h>//
#include <OpenGL/glu.h>//
#include <GLUT/glut.h>
#include<math.h>
#include<stdio.h>
#include <stdlib.h>

#define DEBUG   0
#define MAX_CYLINDERS 100
#define TABLE_SIZE_1 255

/*  Define some constants.  */

#define	PI	3.1415926

/* Max image size allowed. */

#define MAX_SIZE 1024

/*  Define some structures.  */

struct	points	{
    float   x, y, z;
};

typedef struct	rgb_struct	{
    float   r, g, b;
} rgb;

//Noise Table
float noise_tabl[65][65];

/*  Viewing parameters.  */

struct	points	from, at, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz;
float	viewangle, angle, tanv2;
float	xinterval, yinterval;

/*  Illumination parameters.  */

float	lsx, lsy, lsz;
rgb	il, ia;
rgb	ka1, kd1, ks1;
rgb	ka2, kd2, ks2;
rgb	ka3, kd3, ks3;
rgb	ka4, kd4, ks4;
rgb	ka5, kd5, ks5;
rgb	ka6, kd6, ks6;
rgb	tka1, tkd1, tka2, tkd2, tka3, tkd3;
rgb	tka4, tkd4, tka5, tkd5, tka6, tkd6;

int	phong1, phong2, phong3, phong4, phong5, phong6;

/* Cylinder objects. */

int num_cylinders;
int cyl_axis[MAX_CYLINDERS], cyl_phong[MAX_CYLINDERS];
float cyl_x[MAX_CYLINDERS], cyl_y[MAX_CYLINDERS], cyl_z[MAX_CYLINDERS], cyl_r[MAX_CYLINDERS];
rgb cyl_ka[MAX_CYLINDERS], cyl_kd[MAX_CYLINDERS], cyl_ks[MAX_CYLINDERS];
float cyl_min[MAX_CYLINDERS], cyl_max[MAX_CYLINDERS];

/*  Image parameters.  */

int		xmax_pixel, ymax_pixel;

/* Image buffer to stored the ray traced image.  */

float *pixel_rgb;

/*  Object parameters.  */

float	xc1, yc1, zc1, r1, xc2, yc2, zc2, r2, xc3, yc3, zc3, r3;
float	a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
float	px1_min, px1_max, py1_min, py1_max, pz1_min, pz1_max;
float	px2_min, px2_max, py2_min, py2_max, pz2_min, pz2_max;
float	px3_min, px3_max, py3_min, py3_max, pz3_min, pz3_max;

/* Image output file. */

FILE	*outpfile;

/*******************************************************************************
 
 Title:	Read_Information
 
 Purpose:	This function reads in the information about the objects (three spheres
 and three planes) to be rendered.  The information is assumed to be
 stored in an ASCII text file called "ray.dat".
 
 *******************************************************************************/

int Read_Information()
{
    int		i;
    char	str[132], title[132];
    FILE	*inpfile;
    
    if ((inpfile = fopen("khindupur-3.txt", "r")) == NULL) {
        printf("ERROR: Could not open raytrace_demo.txt for read!\n");
        return(0);
    }
    
    
    
    
    
    FILE *meshfile;
    int num_triangles;
    int num_vertices;
    float surf_x_center,surf_y_center,surf_z_center;
    float surf_x_min,surf_y_min,surf_z_min;
    float surf_x_max,surf_y_max,surf_z_max;
    int triangle;
    float tri_normal;
    int vector;
    int vertex;
    float surf_radius;
    int dy,dz;
    
    
    
    
    if ( (meshfile = fopen("meh.txt","rb")) == NULL) {
			     printf("Could not open bsplinemesh.dk for read\n");
        num_triangles = -1;
        num_vertices = -1;
    }
    else {
        
        // Read the triangle mesh data
        
        fread(&surf_x_center, sizeof(float), 1, meshfile);
        fread(&surf_y_center, sizeof(float), 1, meshfile);
        fread(&surf_z_center, sizeof(float), 1, meshfile);
        fread(&surf_x_min, sizeof(float), 1, meshfile);
        fread(&surf_x_max, sizeof(float), 1, meshfile);
        fread(&surf_y_min, sizeof(float), 1, meshfile);
        fread(&surf_y_max, sizeof(float), 1, meshfile);
        fread(&surf_z_min, sizeof(float), 1, meshfile);
        fread(&surf_z_max, sizeof(float), 1, meshfile);
        
        fread(&num_triangles, sizeof(int), 1, meshfile);
        fread(&triangle, sizeof(int), (num_triangles+1)*3, meshfile);
        fread(&tri_normal, sizeof(vector), num_triangles+1, meshfile);
        
        fread(&num_vertices, sizeof(int), 1, meshfile);
        fread(&vertex, sizeof(vector), (num_vertices+1)*3, meshfile);
        fclose(meshfile);
        
        
        printf("Number of triangles in surface mesh: %d\n", num_triangles);
        printf("Number of vertices in surface mesh: %d\n", num_vertices);
        
        printf("Bounding box: (%.2f, %.2f) (%.2f, %.2f) (%.2f, %.2f)\n",
               surf_x_min, surf_x_max, surf_y_min, surf_y_max, surf_z_min, surf_z_max);
        printf("Center: (%.2f %.2f %.2f\n", surf_x_center, surf_y_center, surf_z_center);
        
        // compute the radius of the bounding sphere
        
        surf_radius = surf_x_max - surf_x_min;
        dy = surf_y_max - surf_y_min;
        if (dy > surf_radius) surf_radius = dy;
        dz = surf_z_max - surf_z_min;
        if (dz > surf_radius) surf_radius = dz;
        printf("Bounding sphere radius : %.3f\n", surf_radius);
        
        
        
        
    }
    
    
    
    
    
    /*  Read in viewing information.  */
    
    fscanf(inpfile, "%s %f %f %f\n", str, &from.x, &from.y, &from.z);
    fscanf(inpfile, "%s %f %f %f\n", str, &at.x, &at.y, &at.z);
    fscanf(inpfile, "%s %f %f %f\n", str, &up.x, &up.y, &up.z);
    fscanf(inpfile, "%s %f\n", str, &viewangle);
    angle = viewangle * PI / 180.0;
    tanv2 = tan(angle / 2.0);
    printf("tanv2 = %f\n", tanv2);
    
    /* Read in the viewport information. */
    
    fscanf(inpfile, "%s %f %f\n", str, &VXL, &VXR);
    fscanf(inpfile, "%s %f %f\n", str, &VYB, &VYT);
    
    /* Read in the light vector (L).  Note: To specify a light source position, the light
     vector will need to be computed in Compute_Color().   */
    
    fscanf(inpfile, "%s %f %f %f\n", str, &lsx, &lsy, &lsz);
    
    printf("From: %f %f %f\n", from.x, from.y, from.z);
    printf("At: %f %f %f\n", at.x, at.y, at.z);
    printf("Up: %f %f %f  View Angle: %f \n", up.x, up.y, up.z, viewangle);
    
    /*  Read in spheres' information.  */
    
    fscanf(inpfile, "%s %f %f %f %f\n", str, &xc1, &yc1, &zc1, &r1);
    fscanf(inpfile, "%s %f %f %f %f\n", str, &xc2, &yc2, &zc2, &r2);
    fscanf(inpfile, "%s %f %f %f %f\n", str, &xc3, &yc3, &zc3, &r3);
    
    /*  Read in checker boards' (planes') information.  */
    
    fscanf(inpfile, "%s %f %f %f %f\n", str, &a1, &b1, &c1, &d1);
    fscanf(inpfile, "%s %f %f %f %f\n", str, &a2, &b2, &c2, &d2);
    fscanf(inpfile, "%s %f %f %f %f\n", str, &a3, &b3, &c3, &d3);
    
    /*  Read in the boundaries of the planes.  */
    
    fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px1_min, &py1_min,
           &pz1_min, &px1_max, &py1_max, &pz1_max);
    fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px2_min, &py2_min,
           &pz2_min, &px2_max, &py2_max, &pz2_max);
    fscanf(inpfile, "%s %f %f %f %f %f %f\n", str, &px3_min, &py3_min,
           &pz3_min, &px3_max, &py3_max, &pz3_max);
    
    
    /*  Read in cylinder object information.  */
    
    fscanf(inpfile, "%s %d\n", str, &num_cylinders);
    printf("Number of cylinders: %d\n", num_cylinders);
    if (num_cylinders > MAX_CYLINDERS) {
        printf("Error: num_cylinders exceed MAX_CYLINDERS!\n");
        printf("Program stopped\n");
        return(0);
    }
    
    for (i = 0; i < num_cylinders; i++) {
        fscanf(inpfile, "\n%s\n", title);
        fscanf(inpfile, "%s %d %f %f %f %f\n", str, &cyl_axis[i], &cyl_x[i], &cyl_y[i], &cyl_z[i], &cyl_r[i]);
        fscanf(inpfile, "%s %f %f %f %f %f %f %f %f %f %d\n",
               str, &cyl_ka[i].r, &cyl_ka[i].g, &cyl_ka[i].b,
               &cyl_kd[i].r, &cyl_kd[i].g, &cyl_kd[i].b,
               &cyl_ks[i].r, &cyl_ks[i].g, &cyl_ks[i].b, &cyl_phong[i]);
        fscanf(inpfile, "%s %f %f\n", str, &cyl_min[i], &cyl_max[i]);
        
        printf("Cylinder: %d  Axis %d Center: %.2f %.2f %.2f  Radius: %f\n",
               i, cyl_axis[i], cyl_x[i], cyl_y[i], cyl_z[i], cyl_r[i]);
        printf("            Ka %.1f %.1f %.1f  Kd %.1f %.1f %.1f Ks %.1f %.1f %.1f\n N %d\n",
               cyl_ka[i].r, cyl_ka[i].g, cyl_ka[i].b,
               cyl_kd[i].r, cyl_kd[i].g, cyl_kd[i].b,
               cyl_ks[i].r, cyl_ks[i].g, cyl_ks[i].b, cyl_phong[i]);
        printf("           Min %.2f  Max %.2f \n", cyl_min[i], cyl_max[i]);
    }
    
    
    /*  Read in the image size.  */
    
    fscanf(inpfile, "%s %d %d\n", str, &xmax_pixel, &ymax_pixel);
    
    /* Make sure the image size does not exceed MAX_SIZE x MAX_SIZE.  */
    
    if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE) {
        printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
        printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
        xmax_pixel = MAX_SIZE - 1;
        ymax_pixel = MAX_SIZE - 1;
    }
    
    fclose(inpfile);
    
    /*  Open an output file to store the intensity values of the output image.  */
    
    if ((outpfile = fopen("image.out", "wb")) == NULL) {
        printf("ERROR:  cannot open image.out for write.\n");
        return(0);
    }
    
    /*  Allocate memory for the image buffer.  */
    
    pixel_rgb = new float[3 * xmax_pixel * ymax_pixel];
    printf("image_buf allocated.  Image size %d x %d\n", xmax_pixel, ymax_pixel);
    
    return(1);
}

/*******************************************************************************
 
 Title:	Normalize
 
 Purpose:	This function normalizes the given vector.
 
 *******************************************************************************/

void Normalize(float *x, float *y, float *z)
{
    float	norm;
    
    norm = sqrt(*x * *x + *y * *y + *z * *z);
    if (norm != 0.0) {
        *x = *x / norm;
        *y = *y / norm;
        *z = *z / norm;
    }
}


/*******************************************************************************
 
 Title:	WoodGrain
 
 Purpose:	This function is used to create a woodgrain texture.
 
 *******************************************************************************/

void woodgrain(float u, float v, float w, float *r, float *g, float *b)
{
    
    float radius, angle;
    int grain;
    
    //Creating texture
    radius = sqrt((u*u) + (w*w));
    if (w == 0)
        angle = PI / 2;
    else
        angle = atan2(u, w);
    
    radius = radius + 2 * sin(20 * angle + v / 5);
    grain = (int)radius % 6;
    
    if (grain < 4)
    {
        
        //Setting rgb values
        *r += 0.6*(*r);
        *g += 0.6*(*g);
        *b += 0.6*(*b);
    }
    
    else
    {
        //Setting rgb values
        *r -= 0.6*(*r);
        *g -= 0.6*(*g);
        *b -= 0.6*(*b);
    }
}

/*******************************************************************************
 
 Title:	Calculate Noise
 
 Purpose:	This function is used to create noise for bump map.
 
 *******************************************************************************/

float calc_noise(float iu, float iv, int direction)
{
    int i, j, x, y, left, right;
    float noise, u, v, w, u1, v1, w1;
    float val;
    
    //Filling noise table with random noise
    for (i = 0; i < 65; i++)
        for (j = 0; j < 65; j++)
        {
            
            val = rand() % 255;
            noise_tabl[i][j] = val / 256;
        }
    
    i = (int)iu;
    j = (int)iv;
    x = i%TABLE_SIZE_1;
    y = j%TABLE_SIZE_1;
    
    if (direction == 1)
    {
        if (x < 0)
            left = 0;
        else
            left = x - 1;
        
        if (x >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = x + 1;
        //Getting the value of noise
        noise = (noise_tabl[right][y] - noise_tabl[left][y]) / 2.0;
    }
    
    else {
        if (y <= 0)
            left = 0;
        else
            left = y - 1;
        
        if (y >= TABLE_SIZE_1)
            right = TABLE_SIZE_1;
        else
            right = y + 1;
        //Getting the value of random noise
        noise = (noise_tabl[x][right] - noise_tabl[x][left]) / 2.0;
    }
    return(noise);
    
}


/*******************************************************************************
 
 Title:	Bump Map
 
 Purpose:	This function is used to create bump map texture.
 
 *******************************************************************************/

void Bump_Map(float x, float y, float z, float xc, float yc, float zc, float r, float *nx, float *ny, float *nz)
{
    //Variables
    float xp, yp, zp, iu, iv, xu, yu, zu, xv, yv, zv;
    float fu, fv, a, dx, dy, dz, u, v, nnx, nny, nnz;
    
    xp = (x - xc) / r;
    yp = (y - yc) / r;
    zp = (z - zc) / r;
    
    u = asin(zp);
    v = atan2(yp, xp);
    
    //Getting values to generate the texture
    iu = (u + PI) / (2 * PI) * (TABLE_SIZE_1);
    iv = (v + (PI / 2)) / PI * (TABLE_SIZE_1);
    
    xu = -r*cos(v)*sin(u);
    xv = -r*sin(v)*cos(u);
    yu = r*cos(v)*cos(u);
    yv = -r*sin(v)*sin(u);
    zu = 0.0;
    zv = r*cos(v);
    
    fu = calc_noise(iu, iv, 1);
    fv = calc_noise(iu, iv, 2);
    
    dx = fu*xu + fv*xv;
    dy = fu*yu + fv*yv;
    dz = fu*zu + fv*zv;
    
    
    nnx = *nx;
    nny = *ny;
    nnz = *nz;
    Normalize(&nnx, &nny, &nnz);
    
    //Displacing normals
    dx = fv* (yu*nnz - nny*zu);
    dy = fv* (nnx*zu - xu*nnz);
    dz = fv* (xu*nny - nnx*yu);
    
    dx = fu* (yv*nnz - nny*zv);
    dy = fu* (nnx*zv - xv*nnz);
    dz = fu* (xv*nny - nnx*yv);
    
    Normalize(&dx, &dy, &dz);
    
    a = sqrt(fu*fu + fv*fv);
    dx *= a;
    dy *= a;
    dz *= a;
    
    //Finding new normals after displacement
    *nx += dx;
    *ny += dy;
    *nz += dz;
}



/*******************************************************************************
 
 Title:	Power
 
 Purpose:	This function computes the power of the given base and
 exponent.
 
 *******************************************************************************/

float 	Power(float base, int exp)
{
    int	i;
    float	value;
    
    value = 1.0;
    for (i = 1; i <= exp; i++)
        value *= base;
    
    return(value);
}


/*******************************************************************************
 
 Title:	Compute_M
 
 Purpose:	This function computes the transformation matrix to be used
 in the perspective viewing model.
 
 *******************************************************************************/

void Compute_M()
{
    
    /*  Compute the line-of-sight vector, c.  */
    
    cx = at.x - from.x;
    cy = at.y - from.y;
    cz = at.z - from.z;
    Normalize(&cx, &cy, &cz);
    
    /*  Compute the cross product of vector c and the up vector.  */
    
    ax = cy*up.z - up.y*cz;
    ay = up.x*cz - cx*up.z;
    az = cx*up.y - up.x*cy;
    Normalize(&ax, &ay, &az);
    
    /*  Compute the cross product of vector a and c.  */
    
    bx = ay*cz - cy*az;
    by = cx*az - ax*cz;
    bz = ax*cy - cx*ay;
}

/*******************************************************************************
 
 Title:	Setup_Parameters
 
 Purpose:	This function sets up the necessary parameters for
 performing the ray trace.  It first computes the
 transformation matrix for the perspective viewing model, then
 sets up the default illumination parameters.
 
 *******************************************************************************/

void Setup_Parameters()
{
    
    /*  Compute the transformation matrix for converting world coordinates to eye
     coordinates.  */
    
    Compute_M();
    
    printf("light %f %f %f\n", lsx, lsy, lsz);
    
    /*  Set up the conversion factors for converting from pixel coordinates to
     view port coordinates.  */
    
    xinterval = (VXR - VXL) / xmax_pixel;
    yinterval = (VYT - VYB) / ymax_pixel;
    
    /*  Set up default illumination (Phong lighting) parameters.  */
    
    il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
    ia.r = 1.0;	ia.g = 1.0;	ia.b = 1.0;
    
    /*  Phone lighting parameters for the three spheres.  */
    
    ka1.r = 0.2;	ka1.g = 0.0;	ka1.b = 0.0;
    kd1.r = 0.7;	kd1.g = 0.0;	kd1.b = 0.0;
    ks1.r = 1.0;	ks1.g = 1.0;	ks1.b = 1.0;
    tka1.r = 0.6;	tka1.g = 0.0;	tka1.b = 0.0;
    tkd1.r = 0.6;	tkd1.g = 0.0;	tkd1.b = 0.0;
    phong1 = 60;
    
    ka2.r = 0.0;	ka2.g = 0.2;	ka2.b = 0.0;
    kd2.r = 0.0;	kd2.g = 0.7;	kd2.b = 0.0;
    ks2.r = 1.0;	ks2.g = 1.0;	ks2.b = 1.0;
    tka2.r = 0.0;	tka2.g = 0.6;	tka2.b = 0.0;
    tkd2.r = 0.0;	tkd2.g = 0.6;	tkd2.b = 0.0;
    phong2 = 90;
    
    ka3.r = 0.0;	ka3.g = 0.0;	ka3.b = 0.3;
    kd3.r = 0.0;	kd3.g = 0.0;	kd3.b = 0.7;
    ks3.r = 1.0;	ks3.g = 1.0;	ks3.b = 1.0;
    tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.6;
    tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.6;
    phong3 = 120;
    
    /*  Phone lighting parameters for the three planes.  */
    
    ka4.r = 0.1;	ka4.g = 0.1;	ka4.b = 0.0;
    kd4.r = 0.7;	kd4.g = 0.7;	kd4.b = 0.0;
    ks4.r = 1.0;	ks4.g = 1.0;	ks4.b = 1.0;
    tka4.r = 0.1;	tka4.g = 0.0;	tka4.b = 0.0;
    tkd4.r = 0.7;	tkd4.g = 0.0;	tkd4.b = 0.0;
    phong4 = 120;
    
    ka5.r = 0.1;	ka5.g = 0.0;	ka5.b = 0.1;
    kd5.r = 0.7;	kd5.g = 0.0;	kd5.b = 0.7;
    ks5.r = 1.0;	ks5.g = 1.0;	ks5.b = 1.0;
    tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.1;
    tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;
    phong5 = 120;
    
    ka6.r = 0.1;	ka6.g = 0.1;	ka6.b = 0.1;
    kd6.r = 0.1;	kd6.g = 0.1;	kd6.b = 0.8;
    ks6.r = 1.0;	ks6.g = 1.0;	ks6.b = 1.0;
    tka6.r = 0.0;	tka6.g = 0.1;	tka6.b = 0.0;
    tkd6.r = 0.0;	tkd6.g = 0.7;	tkd6.b = 0.0;
    phong6 = 120;
}


/*******************************************************************************
 
 Title:	Check_Cylinder
 
 Purpose:	This function will check if the give ray intercepts the given
 cylinder.
 
 *******************************************************************************/

void Check_Cylinder(float px, float py, float pz, float dx, float dy, float dz, int axis, float xc, float yc, float zc, float r, float *t1, float *t2)
{
    float	a, b, c, xdiff, ydiff, zdiff, discr;
    
    xdiff = px - xc;
    ydiff = py - yc;
    zdiff = pz - zc;
    
    switch (axis) {
        case 1:
            a = ydiff*ydiff + zdiff*zdiff - r*r;
            b = 2.0*(dy*ydiff + dz*zdiff);
            c = dy*dy + dz*dz;
            break;
            
        case 2:
            a = xdiff*xdiff + zdiff*zdiff - r*r;
            b = 2.0*(dx*xdiff + dz*zdiff);
            c = dx*dx + dz*dz;
            break;
            
        case 3:
            a = xdiff*xdiff + ydiff*ydiff - r*r;
            b = 2.0*(dx*xdiff + dy*ydiff);
            c = dx*dx + dy*dy;
            break;
    }
    
    /*  Check if there are any intersections.  */
    
    discr = b*b - 4.0*a*c;
    if (discr < 0.0) {
        *t1 = -1.0;
        *t2 = -1.0;
    }
    else if (discr == 0.0) {
        *t1 = -b / (2.0*c);
        *t2 = -1.0;
    }
    else {
        discr = sqrt(discr);
        *t1 = (-b + discr) / (2.0*c);
        *t2 = (-b - discr) / (2.0*c);
    }
}



/*******************************************************************************
 
 Title:	Check_Sphere
 
 Purpose:	This function determines if the give ray intercepts the given
 sphere.
 
 *******************************************************************************/

void Check_Sphere(float px, float py, float pz, float dx, float dy, float dz, float xc, float yc, float zc, float r, float *t1, float *t2)
{
    float	a, b, c, xdiff, ydiff, zdiff, discr;
    
    xdiff = px - xc;
    ydiff = py - yc;
    zdiff = pz - zc;
    a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
    b = 2.0*(dx*xdiff + dy*ydiff + dz*zdiff);
    c = dx*dx + dy*dy + dz*dz;
    
    /*  Check if there are any intersections.  */
    
    discr = b*b - 4.0*a*c;
    if (discr < 0.0) {
        *t1 = -1.0;
        *t2 = -1.0;
    }
    else if (discr == 0.0) {
        *t1 = -b / (2.0*c);
        *t2 = -1.0;
    }
    else {
        discr = sqrt(discr);
        *t1 = (-b + discr) / (2.0*c);
        *t2 = (-b - discr) / (2.0*c);
    }
}


/*******************************************************************************
 
 Title:	Check_Plane
 
 Purpose:	This function checks if the given ray intercepts the given
 plane.
 
 *******************************************************************************/

void Check_Plane(float px, float py, float pz, float dx, float dy, float dz, float a, float b, float c, float d, float *t1)
{
    *t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}




/****************************************************************
 Title:check Triangle
 Purpose: this fuction checks if given ray intercepts the triangle
 
 
 **************************************************************/


float dot(float v0[3], float v1[3])
{
    return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] );
}

#define MAX_VERTICES 2000
#define MAX_TRIANGLES 1000

float vertex[MAX_VERTICES][3];
int triangle[MAX_TRIANGLES][3];

int tri;

void Check_Triangle(int triNumber, float px,float py,float pz,float dx,float dy,float dz,float a,float b,float c,float d,float *t1)
{
    int	 vert1, vert2, vert3;
    float dot00, dot01, dot02, dot11, dot12, invDenom;
    float t, ipx, ipy, ipz, u, v, p[3][3], v0[3], v1[3], v2[3];
    
    t = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
    if (t > 0.00001) {
        
        // Check if the intersection point is inside the triangle.
        
        ipx = px + t*dx;
        ipy = py + t*dy;
        ipz = pz + t*dz;
        
        vert1 = triangle[tri][0];
        vert2 = triangle[tri][1];
        vert3 = triangle[tri][2];
        
        p[0][0] = vertex[vert1][0];
        p[0][1] = vertex[vert1][1];
        p[0][2] = vertex[vert1][2];
        
        p[1][0] = vertex[vert2][0];
        p[1][1] = vertex[vert2][1];
        p[1][2] = vertex[vert2][2];
        
        p[2][0] = vertex[vert3][0];
        p[2][1] = vertex[vert3][1];
        p[2][2] = vertex[vert3][2];
        
        v0[0] = p[2][0] - p[0][0];
        v0[1] = p[2][1] - p[0][1];
        v0[2] = p[2][2] - p[0][2];
        
        v1[0] = p[1][0] - p[0][0];
        v1[1] = p[1][1] - p[0][1];
        v1[2] = p[1][2] - p[0][2];
        
        v2[0] = ipx - p[0][0];
        v2[1] = ipy - p[0][1];
        v2[2] = ipz - p[0][2];
        
        dot00 = dot(v0, v0);
        dot01 = dot(v0, v1);
        dot02 = dot(v0, v2);
        dot11 = dot(v1, v1);
        dot12 = dot(v1, v2);
        
        invDenom = 1.0/( dot00*dot11 - dot01*dot01);
        
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;
        
        if ( u>=0.0 && v >= 0.0 && u+v <= 1.0)
            *t1 = t;
        else
            *t1 = -1;
    }
    else
        *t1 = -1;
}



/*******************************************************************************
 
 Title:	Compute_Intersection
 
 Purpose:	This function computes the intersection of ray with an
 object.  The intersection point is given by a parametric value
 t, where ray = p + d*t, d = the direction of the ray, and p is
 the starting point of the ray.
 
 *******************************************************************************/

void Compute_Intersection(float px, float py, float pz, float dx, float dy, float dz, float t, float *newx, float *newy, float *newz)
{
    *newx = px + t*dx;
    *newy = py + t*dy;
    *newz = pz + t*dz;
}


/*******************************************************************************
 
 Title:	Compute_Color
 
 Purpose:	This function computes the intensity of the color for the
 given location based on the Phong lighting model.
 
 *******************************************************************************/

void Compute_Color(int shadow_flag, float ipx, float ipy, float  ipz, float  nx, float  ny, float  nz,
                   rgb ia, rgb ka, rgb kd, rgb ks, int n, float *r, float *g, float *b, int obj_num)
{
    float	vx, vy, vz, rx, ry, rz;
    float	ndotl, vdotr, cosalphapower;
    float lx, ly, lz;
    /*  Compute the view vector.  */
    
    vx = from.x - ipx;
    vy = from.y - ipy;
    vz = from.z - ipz;
    Normalize(&vx, &vy, &vz);
    
    lx = lsx - ipx;
    ly = lsy - ipy;
    lz = lsz - ipz;
    Normalize(&lx, &ly, &lz);
    
    /*  Compute the R (reflection) vector.  */
    
    ndotl = nx*lx + ny*ly + nz*lz;
    rx = -(2.0*ndotl*nx - lx);
    ry = -(2.0*ndotl*ny - ly);
    rz = -(2.0*ndotl*nz - lz);
    
    /* Compute the V (view) vector. */
    
    vdotr = vx*rx + vy*ry + vz*rz;
    
    //If plane 3 , apply woodgrain pattern
    if (obj_num == 6)
    {
        *g = *g + tka5.g*ndotl*il.g;
        *b = *b + tkd5.b*ndotl*il.b;
        woodgrain(ipx, ipy, ipz, r, g, b);
    }
    
    /* Compute Ia * Ka.  */
    else
    {
        
        *r = ia.r * ka.r;
        *g = ia.g * ka.g;
        *b = ia.b * ka.b;
        
        /* Compute diffuse reflection. */
        
        if (ndotl >= 0.0 && shadow_flag == 0) {
            
            /*  diffuse reflection = kd * N dot L * Il  */
            
            *r = *r + kd.r*ndotl*il.r;
            *g = *g + kd.g*ndotl*il.g;
            *b = *b + kd.b*ndotl*il.b;
            
            if (vdotr >= 0.0) {
                
                /*  specular reflection = ks * cos(alpha)**K^n * Il */
                
                cosalphapower = Power(vdotr, n);
                
                *r = *r + ks.r*cosalphapower*il.r;
                *g = *g + ks.g*cosalphapower*il.g;
                *b = *b + ks.b*cosalphapower*il.b;
            }
        }
    }
    if (shadow_flag == 1)
    {
        //If there is a shadow, make that rgb black
        *r = 0.0;
        *g = 0.0;
        *b = 0.0;
        
    }
    
    /*  Make sure that the color is within range.  */
    
    if (*r > 1.0) *r = 1.0;
    if (*g > 1.0) *g = 1.0;
    if (*b > 1.0) *b = 1.0;
}


/*******************************************************************************
 
 Title:	Reflection
 
 Purpose:	This function computes the reflected ray and returns rgb values.
 
 *******************************************************************************/

void Reflection(float fromx, float fromy, float fromz, float tox, float toy, float toz, float *rr, float *rg, float *rb, int object_id)
{
    int	    xp, yp, obj_num, shadow_flag = 0;
    int	    i, texture, buf_ptr, index, cyl_id;
    int     num_image = 1;
    float	xv, yv, dx, dy, dz, nx, ny, nz;
    float	t_min, t1, t2;
    float ipx, ipy, ipz;
    float	r, g, b;
    float   u, v;
    
    
    t_min = 999.0;
    obj_num = 0;
    texture = 0;
    
    /*  Check if the current ray intercepts the cyinders.  */
    
    for (i = 0; i<num_cylinders; i++) {
        if (object_id == 7 + i) continue;
        Check_Cylinder(fromx, fromy, fromz, tox, toy, toz, cyl_axis[i], cyl_x[i], cyl_y[i], cyl_z[i],
                       cyl_r[i], &t1, &t2);
        
        if (t1>0.00001 && t1 < t_min) {
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i])
                        t1 = -1;
                    break;
            }
            
            if (t1 > 0.0) {
                t_min = t1;
                cyl_id = i;
                
                // cylinder objects' number start at #7 since object number #6 is plane #3.
                obj_num = 7 + i;
                
                // compute the surface normal
                switch (cyl_axis[i]) {
                    case 1:
                        nx = 0.0;
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 2:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = 0.0;
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 3:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = 0.0;
                        break;
                }
            }
        }
        
        if (t2 > 0.00001 && t2<t_min) {
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
            }
            
            if (t2 > 0.0) {
                t_min = t2;
                // cylinder objects' number start at #7 since object number #6 is plane #3
                obj_num = 7 + i;
                cyl_id = i;
                
                // compute the surface normal
                switch (cyl_axis[i]) {
                    case 1:
                        nx = 0.0;
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 2:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = 0.0;
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 3:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = 0.0;
                        break;
                }
            }
        }
    }
    
    /*  Check if the current ray intercepts sphere #1.  */
    
    if (object_id != 1){
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc1,
                     yc1, zc1, r1, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 1;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 1;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    /*  Check if the current ray intercepts sphere #2.  */
    if (object_id != 2){
        
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc2,
                     yc2, zc2, r2, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 2;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 2;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    /*  Check if the current ray intercepts sphere #3.  */
    
    if (object_id != 3){
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc3,
                     yc3, zc3, r3, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 3;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 3;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    
    /*  Check if the current ray intercepts plane #1.  */
    if (object_id != 4){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a1, b1, c1, d1, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px1_min && ipx <= px1_max &&
                ipy >= py1_min  && ipy <= py1_max &&
                ipz >= pz1_min && ipz <= pz1_max) {
                
                t_min = t1;
                obj_num = 4;
            }
        }
        
    }
    /*  Check if the current ray intercepts plane #2.  */
    
    if (object_id != 5){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a2, b2, c2, d2, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px2_min && ipx <= px2_max &&
                ipy >= py2_min  && ipy <= py2_max &&
                ipz >= pz2_min && ipz <= pz2_max) {
                
                t_min = t1;
                obj_num = 5;
            }
        }
        
    }
    /*  Check if the current ray intercepts plane #3.  */
    
    if (object_id != 6){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a3, b3, c3, d3, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px3_min && ipx <= px3_max &&
                ipy >= py3_min  && ipy <= py3_max &&
                ipz >= pz3_min && ipz <= pz3_max) {
                
                t_min = t1;
                obj_num = 6;
            }
        }
        
        
    }
    /*  Compute the intensity to use at the current pixel.  */
    
    switch (obj_num) {
            
            /*  The current ray does not intersect any of the objects.  */
            
        case 0: r = 0.0;
            g = 0.0;
            b = 0.0;
            break;
            
            /*  The current ray intercept sphere #1.  */
            
        case 1:
            nx = ipx - xc1;
            ny = ipy - yc1;
            nz = ipz - zc1;
            Normalize(&nx, &ny, &nz);
            
            
            texture = 0;
            
            Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b, object_id);
            break;
            
            /*  The current ray intercepts sphere #2.  */
            
        case 2:
            nx = ipx - xc2;
            ny = ipy - yc2;
            nz = ipz - zc2;
            Normalize(&nx, &ny, &nz);
            
            texture = 0;
            Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts sphere #3.  */
            
        case 3:
            nx = ipx - xc3;
            ny = ipy - yc3;
            nz = ipz - zc3;
            Normalize(&nx, &ny, &nz);
            
            texture = 0;
            
            // Compute texture. */
            
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks3, phong3, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #1.  */
            
        case 4:
            nx = a1;
            ny = b1;
            nz = c1;
            
            
            if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4, tkd4, ks4, phong4, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #2.  */
            
        case 5:
            nx = a2;
            ny = b2;
            nz = c2;
            
            
            if (ipz < 2.0 || (ipz >= 4.0 && ipz < 6.0)) {
                if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong5, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka5, kd5, ks5, phong5, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #3.  */
            
        case 6:
            nx = a3;
            ny = b3;
            nz = c3;
            
            
            if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                if ((ipz >= 2.0 && ipz < 4.0) || (ipz >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipz < 2.0) || (ipz >= 4.0 && ipz < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b, obj_num);
            break;
    }
    
    // Compute the color if the intersected object is a cylinder.
    if (obj_num > 6) {
        
        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, cyl_ka[cyl_id], cyl_kd[cyl_id], cyl_ks[cyl_id],
                      cyl_phong[cyl_id], &r, &g, &b, obj_num);
    }
    //Getting new rgb values
    *rr = r;
    *rg = g;
    *rb = b;
    
}

/*******************************************************************************
 
 Title:	Refraction
 
 Purpose:	This function computes the refraction and returns rgb values.
 
 
 
 *******************************************************************************/


void Refraction(float fromx, float fromy, float fromz, float tox, float toy, float toz, float *rr, float *rg, float *rb, int object_id)
{
    int	    xp, yp, obj_num, shadow_flag = 0;
    int	    i, texture, buf_ptr, index, cyl_id;
    int     num_image = 1;
    float	xv, yv, dx, dy, dz, nx, ny, nz;
    float	t_min, t1, t2;
    float ipx, ipy, ipz;
    float	r, g, b;
    float   u, v;
    
    
    t_min = 999.0;
    obj_num = 0;
    texture = 0;
    
    /*  Check if the current ray intercepts the cyinders.  */
    
    for (i = 0; i<num_cylinders; i++) {
        if (object_id == 7 + i) continue;
        Check_Cylinder(fromx, fromy, fromz, tox, toy, toz, cyl_axis[i], cyl_x[i], cyl_y[i], cyl_z[i],
                       cyl_r[i], &t1, &t2);
        
        if (t1>0.00001 && t1 < t_min) {
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i])
                        t1 = -1;
                    break;
            }
            
            if (t1 > 0.0) {
                t_min = t1;
                cyl_id = i;
                
                // cylinder objects' number start at #7 since object number #6 is plane #3.
                obj_num = 7 + i;
                
                // compute the surface normal
                switch (cyl_axis[i]) {
                    case 1:
                        nx = 0.0;
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 2:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = 0.0;
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 3:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = 0.0;
                        break;
                }
            }
        }
        
        if (t2 > 0.00001 && t2<t_min) {
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
            }
            
            if (t2 > 0.0) {
                t_min = t2;
                // cylinder objects' number start at #7 since object number #6 is plane #3
                obj_num = 7 + i;
                cyl_id = i;
                
                // compute the surface normal
                switch (cyl_axis[i]) {
                    case 1:
                        nx = 0.0;
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 2:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = 0.0;
                        nz = (ipz - cyl_z[i]) / cyl_r[i];
                        break;
                        
                    case 3:
                        nx = (ipx - cyl_x[i]) / cyl_r[i];
                        ny = (ipy - cyl_y[i]) / cyl_r[i];
                        nz = 0.0;
                        break;
                }
            }
        }
    }
    
    /*  Check if the current ray intercepts sphere #1.  */
    
    if (object_id != 1){
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc1,
                     yc1, zc1, r1, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 1;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 1;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    /*  Check if the current ray intercepts sphere #2.  */
    if (object_id != 2){
        
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc2,
                     yc2, zc2, r2, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 2;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 2;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    /*  Check if the current ray intercepts sphere #3.  */
    
    if (object_id != 3){
        Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc3,
                     yc3, zc3, r3, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            t_min = t1;
            obj_num = 3;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            t_min = t2;
            obj_num = 3;
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t2, &ipx, &ipy, &ipz);
        }
        
    }
    
    /*  Check if the current ray intercepts plane #1.  */
    if (object_id != 4){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a1, b1, c1, d1, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px1_min && ipx <= px1_max &&
                ipy >= py1_min  && ipy <= py1_max &&
                ipz >= pz1_min && ipz <= pz1_max) {
                
                t_min = t1;
                obj_num = 4;
            }
        }
        
    }
    /*  Check if the current ray intercepts plane #2.  */
    
    if (object_id != 5){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a2, b2, c2, d2, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px2_min && ipx <= px2_max &&
                ipy >= py2_min  && ipy <= py2_max &&
                ipz >= pz2_min && ipz <= pz2_max) {
                
                t_min = t1;
                obj_num = 5;
            }
        }
        
    }
    /*  Check if the current ray intercepts plane #3.  */
    
    if (object_id != 6){
        Check_Plane(fromx, fromy, fromz, tox, toy, toz, a3, b3, c3, d3, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(fromx, fromy, fromz,
                                 tox, toy, toz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px3_min && ipx <= px3_max &&
                ipy >= py3_min  && ipy <= py3_max &&
                ipz >= pz3_min && ipz <= pz3_max) {
                
                t_min = t1;
                obj_num = 6;
            }
        }
        
        
    }
    /*  Compute the intensity to use at the current pixel.  */
    
    switch (obj_num) {
            
            /*  The current ray does not intersect any of the objects.  */
            
        case 0: r = 0.0;
            g = 0.0;
            b = 0.0;
            break;
            
            /*  The current ray intercept sphere #1.  */
            
        case 1:
            nx = ipx - xc1;
            ny = ipy - yc1;
            nz = ipz - zc1;
            Normalize(&nx, &ny, &nz);
            
            
            texture = 0;
            
            Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b, object_id);
            break;
            
            /*  The current ray intercepts sphere #2.  */
            
        case 2:
            nx = ipx - xc2;
            ny = ipy - yc2;
            nz = ipz - zc2;
            Normalize(&nx, &ny, &nz);
            
            texture = 0;
            Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts sphere #3.  */
            
        case 3:
            nx = ipx - xc3;
            ny = ipy - yc3;
            nz = ipz - zc3;
            Normalize(&nx, &ny, &nz);
            
            texture = 0;
            
            // Compute texture. */
            
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks3, phong3, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #1.  */
            
        case 4:
            nx = a1;
            ny = b1;
            nz = c1;
            
            
            if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4, tkd4, ks4, phong4, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #2.  */
            
        case 5:
            nx = a2;
            ny = b2;
            nz = c2;
            
            
            if (ipz < 2.0 || (ipz >= 4.0 && ipz < 6.0)) {
                if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong5, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka5, kd5, ks5, phong5, &r, &g, &b, obj_num);
            break;
            
            /*  The current ray intercepts checker board #3.  */
            
        case 6:
            nx = a3;
            ny = b3;
            nz = c3;
            
            
            if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                if ((ipz >= 2.0 && ipz < 4.0) || (ipz >= 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            else {
                if ((ipz < 2.0) || (ipz >= 4.0 && ipz < 6.0))
                    texture = 1;
                else
                    texture = 0;
            }
            
            if (texture == 1)
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &r, &g, &b, obj_num);
            else
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b, obj_num);
            break;
    }
    
    // Compute the color if the intersected object is a cylinder.
    if (obj_num > 6) {
        
        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, cyl_ka[cyl_id], cyl_kd[cyl_id], cyl_ks[cyl_id],
                      cyl_phong[cyl_id], &r, &g, &b, obj_num);
    }
    
    //Setting new rgb values
    *rr = r;
    *rg = g;
    *rb = b;
    
}

/*******************************************************************************
 
 Title:	Check Shadow
 
 Purpose:	This function checks for intersection to see if a shadow is generated
 and returns 1 or 0 (shadow flag value)
 
 *******************************************************************************/


int Check_Shadow(float x, float y, float z, int object_id)
{
    float t1, t2, t_min, ipx, ipy, ipz;
    int i;
    
    t_min = 999.0;
    
    /*  Check if the current ray intercepts the cyinders.  */
    
    for (i = 0; i<num_cylinders; i++) {
        if (object_id == 7 + i) continue;
        Check_Cylinder(x, y, z, lsx, lsy, lsz, cyl_axis[i], cyl_x[i], cyl_y[i], cyl_z[i],
                       cyl_r[i], &t1, &t2);
        
        if (t1>0.00001 && t1 < t_min) {
            
            Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i])
                        t1 = -1;
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i])
                        t1 = -1;
                    break;
            }
            
            if (t1 > 0.0) {
                return 1;
            }
        }
        
        if (t2 > 0.00001 && t2<t_min) {
            Compute_Intersection(x, y, z, lsx, lsy, lsz, t2, &ipx, &ipy, &ipz);
            
            switch (cyl_axis[i]) {
                case 1:
                    if (ipx < cyl_min[i] || ipx > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 2:
                    if (ipy < cyl_min[i] || ipy > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
                    
                case 3:
                    if (ipz < cyl_min[i] || ipz > cyl_max[i]) {
                        if (t1 > 0.0)
                            t2 = t1;
                        else
                            t2 = -1;
                    }
                    break;
            }
            
            if (t2 > 0.0) {
                return 1;
            }
        }
    }
    
    /*  Check if the current ray intercepts sphere #1.  */
    
    if (object_id != 1){
        Check_Sphere(x, y, z, lsx, lsy, lsz, xc1,
                     yc1, zc1, r1, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            return 1;
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            return 1;
        }
    }
    /*  Check if the current ray intercepts sphere #2.  */
    if (object_id != 2){
        Check_Sphere(x, y, z, lsx, lsy, lsz, xc2,
                     yc2, zc2, r2, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            return 1;
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            return 1;
        }
    }
    /*  Check if the current ray intercepts sphere #3.  */
    if (object_id != 3){
        Check_Sphere(x, y, z, lsx, lsy, lsz, xc3,
                     yc3, zc3, r3, &t1, &t2);
        
        if (t1 >= 0.0 && t1 < t_min) {
            return 1;
        }
        
        if (t2 >= 0.0 && t2 < t_min) {
            return 1;
        }
    }
    /*  Check if the current ray intercepts plane #1.  */
    if (object_id != 4){
        Check_Plane(x, y, z, lsx, lsy, lsz, a1, b1, c1, d1, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px1_min && ipx <= px1_max &&
                ipy >= py1_min  && ipy <= py1_max &&
                ipz >= pz1_min && ipz <= pz1_max) {
                return 1;
            }
        }
    }
    /*  Check if the current ray intercepts plane #2.  */
    if (object_id != 5){
        Check_Plane(x, y, z, lsx, lsy, lsz, a2, b2, c2, d2, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px2_min && ipx <= px2_max &&
                ipy >= py2_min  && ipy <= py2_max &&
                ipz >= pz2_min && ipz <= pz2_max) {
                return 1;
            }
        }
    }
    /*  Check if the current ray intercepts plane #3.  */
    if (object_id != 6){
        Check_Plane(x, y, z, lsx, lsy, lsz, a3, b3, c3, d3, &t1);
        
        if (t1 >= 0.0 && t1 < t_min) {
            /*  Check if the intersection point is inside the min/max values. */
            
            Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);
            
            if (ipx >= px3_min && ipx <= px3_max &&
                ipy >= py3_min  && ipy <= py3_max &&
                ipz >= pz3_min && ipz <= pz3_max) {
                return 1;
            }
        }
    }
    return 0;
}


/*******************************************************************************
 
 Title:	Ray_Trace
 
 Purpose:	This function performs simple ray tracing.  *******************************************************************************/

void Ray_Trace()
{
    int	    xp, yp, obj_num, shadow_flag;
    int	    i, texture, buf_ptr, index, cyl_id;
    int     num_image = 1;
    float	xv, yv, dx, dy, dz, nx, ny, nz;
    float	t_min, t1, t2, ipx, ipy, ipz;
    float	r, g, b;
    float   rr, rg, rb;
    float   rrr, rgg, rbb;
    float   interx, intery, interz;
    float   reflectedx, reflectedy, reflectedz;
    float   dot;
    float   u, v;
    float   n1 = 1.0;
    float   n2 = 1.4;
    float   transx, transy, transz;
    
    
    /*  Generate a ray for each pixel in the desired image.  */
    
    printf("Ray tracing...\n");
    buf_ptr = 0;
    for (xp = 0; xp < xmax_pixel; xp++) {
        u = (float)xp / xmax_pixel;
        
        for (yp = 0; yp < ymax_pixel; yp++) {
            v = (float)yp / ymax_pixel;
            
            /*  Compute the corresponding view port coordinates.  */
            
            xv = VXL + xp * xinterval;
            yv = VYB + yp * yinterval;
            
            /*  Compute the direction of the current ray from the "From" point to the
             current position on the image.  */
            
            dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
            dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
            dz = az*xv*tanv2 + bz*yv*tanv2 + cz;
            
            
            t_min = 999.0;
            obj_num = 0;
            texture = 0;
            
            /*  Check if the current ray intercepts the cyinders.  */
            
            for (i = 0; i<num_cylinders; i++) {
                
                Check_Cylinder(from.x, from.y, from.z, dx, dy, dz, cyl_axis[i], cyl_x[i], cyl_y[i], cyl_z[i],
                               cyl_r[i], &t1, &t2);
                
                if (t1>0.00001 && t1 < t_min) {
                    
                    Compute_Intersection(from.x, from.y, from.z,
                                         dx, dy, dz, t1, &ipx, &ipy, &ipz);
                    
                    switch (cyl_axis[i]) {
                        case 1:
                            if (ipx < cyl_min[i] || ipx > cyl_max[i])
                                t1 = -1;
                            break;
                            
                        case 2:
                            if (ipy < cyl_min[i] || ipy > cyl_max[i])
                                t1 = -1;
                            break;
                            
                        case 3:
                            if (ipz < cyl_min[i] || ipz > cyl_max[i])
                                t1 = -1;
                            break;
                    }
                    
                    if (t1 > 0.0) {
                        t_min = t1;
                        cyl_id = i;
                        
                        // cylinder objects' number start at #7 since object number #6 is plane #3.
                        obj_num = 7 + i;
                        
                        // compute the surface normal
                        switch (cyl_axis[i]) {
                            case 1:
                                nx = 0.0;
                                ny = (ipy - cyl_y[i]) / cyl_r[i];
                                nz = (ipz - cyl_z[i]) / cyl_r[i];
                                break;
                                
                            case 2:
                                nx = (ipx - cyl_x[i]) / cyl_r[i];
                                ny = 0.0;
                                nz = (ipz - cyl_z[i]) / cyl_r[i];
                                break;
                                
                            case 3:
                                nx = (ipx - cyl_x[i]) / cyl_r[i];
                                ny = (ipy - cyl_y[i]) / cyl_r[i];
                                nz = 0.0;
                                break;
                        }
                    }
                }
                
                if (t2 > 0.00001 && t2<t_min) {
                    Compute_Intersection(from.x, from.y, from.z,
                                         dx, dy, dz, t2, &ipx, &ipy, &ipz);
                    
                    switch (cyl_axis[i]) {
                        case 1:
                            if (ipx < cyl_min[i] || ipx > cyl_max[i]) {
                                if (t1 > 0.0)
                                    t2 = t1;
                                else
                                    t2 = -1;
                            }
                            break;
                            
                        case 2:
                            if (ipy < cyl_min[i] || ipy > cyl_max[i]) {
                                if (t1 > 0.0)
                                    t2 = t1;
                                else
                                    t2 = -1;
                            }
                            break;
                            
                        case 3:
                            if (ipz < cyl_min[i] || ipz > cyl_max[i]) {
                                if (t1 > 0.0)
                                    t2 = t1;
                                else
                                    t2 = -1;
                            }
                            break;
                    }
                    
                    if (t2 > 0.0) {
                        t_min = t2;
                        // cylinder objects' number start at #7 since object number #6 is plane #3
                        obj_num = 7 + i;
                        cyl_id = i;
                        
                        // compute the surface normal
                        switch (cyl_axis[i]) {
                            case 1:
                                nx = 0.0;
                                ny = (ipy - cyl_y[i]) / cyl_r[i];
                                nz = (ipz - cyl_z[i]) / cyl_r[i];
                                break;
                                
                            case 2:
                                nx = (ipx - cyl_x[i]) / cyl_r[i];
                                ny = 0.0;
                                nz = (ipz - cyl_z[i]) / cyl_r[i];
                                break;
                                
                            case 3:
                                nx = (ipx - cyl_x[i]) / cyl_r[i];
                                ny = (ipy - cyl_y[i]) / cyl_r[i];
                                nz = 0.0;
                                break;
                        }
                    }
                }
            }
            
            /*  Check if the current ray intercepts sphere #1.  */
            
            
            Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc1,
                         yc1, zc1, r1, &t1, &t2);
            
            if (t1 >= 0.0 && t1 < t_min) {
                t_min = t1;
                obj_num = 1;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
            }
            
            if (t2 >= 0.0 && t2 < t_min) {
                t_min = t2;
                obj_num = 1;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t2, &ipx, &ipy, &ipz);
            }
            
            /*  Check if the current ray intercepts sphere #2.  */
            
            Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc2,
                         yc2, zc2, r2, &t1, &t2);
            
            if (t1 >= 0.0 && t1 < t_min) {
                t_min = t1;
                obj_num = 2;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
            }
            
            if (t2 >= 0.0 && t2 < t_min) {
                t_min = t2;
                obj_num = 2;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t2, &ipx, &ipy, &ipz);
            }
            
            /*  Check if the current ray intercepts sphere #3.  */
            
            Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc3,
                         yc3, zc3, r3, &t1, &t2);
            
            if (t1 >= 0.0 && t1 < t_min) {
                t_min = t1;
                obj_num = 3;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
            }
            
            if (t2 >= 0.0 && t2 < t_min) {
                t_min = t2;
                obj_num = 3;
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t2, &ipx, &ipy, &ipz);
            }
            
            /*  Check if the current ray intercepts plane #1.  */
            
            Check_Plane(from.x, from.y, from.z, dx, dy, dz, a1, b1, c1, d1, &t1);
            
            if (t1 >= 0.0 && t1 < t_min) {
                /*  Check if the intersection point is inside the min/max values. */
                
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
                
                if (ipx >= px1_min && ipx <= px1_max &&
                    ipy >= py1_min  && ipy <= py1_max &&
                    ipz >= pz1_min && ipz <= pz1_max) {
                    
                    t_min = t1;
                    obj_num = 4;
                }
            }
            
            /*  Check if the current ray intercepts plane #2.  */
            
            Check_Plane(from.x, from.y, from.z, dx, dy, dz, a2, b2, c2, d2, &t1);
            
            if (t1 >= 0.0 && t1 < t_min) {
                /*  Check if the intersection point is inside the min/max values. */
                
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
                
                if (ipx >= px2_min && ipx <= px2_max &&
                    ipy >= py2_min  && ipy <= py2_max &&
                    ipz >= pz2_min && ipz <= pz2_max) {
                    
                    t_min = t1;
                    obj_num = 5;
                }
            }
            
            /*  Check if the current ray intercepts plane #3.  */
            
            Check_Plane(from.x, from.y, from.z, dx, dy, dz, a3, b3, c3, d3, &t1);
            
            if (t1 >= 0.0 && t1 < t_min) {
                /*  Check if the intersection point is inside the min/max values. */
                
                Compute_Intersection(from.x, from.y, from.z,
                                     dx, dy, dz, t1, &ipx, &ipy, &ipz);
                
                if (ipx >= px3_min && ipx <= px3_max &&
                    ipy >= py3_min  && ipy <= py3_max &&
                    ipz >= pz3_min && ipz <= pz3_max) {
                    
                    t_min = t1;
                    obj_num = 6;
                }
            }
            
            
            
            /*  Compute the intensity to use at the current pixel.  */
            
            switch (obj_num) {
                    
                    /*  The current ray does not intersect any of the objects.  */
                    
                case 0: r = 0.0;
                    g = 0.7;
                    b = 0.8;
                    
                    rr = 0.0;
                    rg = 0.7;
                    rb = 0.8;
                    
                    rrr = 0.0;
                    rgg = 0.7;
                    rbb = 0.8;
                    
                    break;
                    
                    /*  The current ray intercept sphere #1.  */
                    
                case 1:
                    nx = ipx - xc1;
                    ny = ipy - yc1;
                    nz = ipz - zc1;
                    Normalize(&nx, &ny, &nz);
                    
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    texture = 0;
                    
                    //Finding intersection ray
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    //Finding dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b, obj_num);
                    break;
                    
                    /*  The current ray intercepts sphere #2.  */
                    
                case 2:
                    nx = ipx - xc2;
                    ny = ipy - yc2;
                    nz = ipz - zc2;
                    Normalize(&nx, &ny, &nz);
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    
                    //Finding intersection ray
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    //Finding dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    texture = 0;
                    
                    Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b, obj_num);
                    break;
                    
                    /*  The current ray intercepts sphere #3.  */
                    
                case 3:
                    nx = ipx - xc3;
                    ny = ipy - yc3;
                    nz = ipz - zc3;
                    Normalize(&nx, &ny, &nz);
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    
                    //Finding intersection ray
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    //Finding dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    // Compute texture. */
                    
                    if (texture == 1){
                        Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz);
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks3, phong3, &r, &g, &b, obj_num);
                    }
                    
                    else{
                        Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz);
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &r, &g, &b, obj_num);
                    }
                    break;
                    
                    /*  The current ray intercepts checker board #1.  */
                    
                case 4:
                    nx = a1;
                    ny = b1;
                    nz = c1;
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    
                    //Finding intersection ray
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    
                    //Finding dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                        if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    else {
                        if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    if (texture == 1){
                        
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4, tkd4, ks4, phong4, &r, &g, &b, obj_num);
                    }
                    else{
                        
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b, obj_num);
                    }
                    break;
                    
                    /*  The current ray intercepts checker board #2.  */
                    
                case 5:
                    nx = a2;
                    ny = b2;
                    nz = c2;
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    
                    //Finding intersection of ray
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    
                    //Finding dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    
                    if (ipz < 2.0 || (ipz >= 4.0 && ipz < 6.0)) {
                        if ((ipy >= 2.0 && ipy < 4.0) || (ipy >= 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    else {
                        if ((ipy < 2.0) || (ipy >= 4.0 && ipy < 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    if (texture == 1){
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong5, &r, &g, &b, obj_num);
                    }
                    else{
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka5, kd5, ks5, phong5, &r, &g, &b, obj_num);
                    }
                    break;
                    
                    /*  The current ray intercepts checker board #3.  */
                    
                case 6:
                    nx = a3;
                    ny = b3;
                    nz = c3;
                    
                    shadow_flag = 0;
                    shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                    
                    //Finding intersection point
                    interx = dx - ipx;
                    intery = dy - ipy;
                    interz = dz - ipz;
                    
                    //Computing dot product and reflected ray
                    dot = nx*interx + ny*intery + nz*interz;
                    reflectedx = interx - 2 * (nx*dot);
                    reflectedy = intery - 2 * (ny*dot);
                    reflectedz = interz - 2 * (nz*dot);
                    
                    Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                    
                    //Formula for refraction
                    transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                    transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                    transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                    
                    Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                    
                    
                    if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
                        if ((ipz >= 2.0 && ipz < 4.0) || (ipz >= 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    else {
                        if ((ipz < 2.0) || (ipz >= 4.0 && ipz < 6.0))
                            texture = 1;
                        else
                            texture = 0;
                    }
                    
                    if (texture == 1)
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &r, &g, &b, obj_num);
                    else
                        Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b, obj_num);
                    break;
            }
            
            // Compute the color if the intersected object is a cylinder.
            if (obj_num > 6) {
                shadow_flag = 0;
                
                shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
                
                interx = dx - ipx;
                intery = dy - ipy;
                interz = dz - ipz;
                
                dot = nx*interx + ny*intery + nz*interz;
                reflectedx = interx - 2 * (nx*dot);
                reflectedy = intery - 2 * (ny*dot);
                reflectedz = interz - 2 * (nz*dot);
                
                Reflection(ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &rr, &rg, &rb, obj_num);
                
                transx = ((n1 / n2)*interx) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nx;
                transy = ((n1 / n2)*intery) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*ny;
                transz = ((n1 / n2)*interz) - ((sqrt(1 - ((n1*n1) / ((n2*n2))*(1 - (dot*dot))))) + (n1 / n2)*dot)*nz;
                
                Refraction(ipx, ipy, ipz, transx, transy, transz, &rrr, &rgg, &rbb, obj_num);
                
                Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, cyl_ka[cyl_id], cyl_kd[cyl_id], cyl_ks[cyl_id],
                              cyl_phong[cyl_id], &r, &g, &b, obj_num);
            }
            /* Save the computed color intensity to the image buffer. */
            index = 3 * (xp + xmax_pixel * yp);
            if (obj_num > 6){
                
                pixel_rgb[index] = r + 0.3*rr;
                pixel_rgb[index + 1] = g + 0.3*rg;
                pixel_rgb[index + 2] = b + 0.3*rb;
            }
            else if (obj_num == 1 || obj_num == 2)
            {
                pixel_rgb[index] = r + 0.8*rrr;
                pixel_rgb[index + 1] = g + 0.8*rgg;
                pixel_rgb[index + 2] = b + 0.8*rbb;
            }
            else
            {
                
                pixel_rgb[index] = r;
                pixel_rgb[index + 1] = g;
                pixel_rgb[index + 2] = b;
            }
        }
    }
    
    /*  Write the image to the output file.  */
    
    printf("Writing to image...\n");
    fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
    fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
    fwrite(pixel_rgb, sizeof(float), 3 * xmax_pixel*ymax_pixel, outpfile);
    fclose(outpfile);
}


/* Initialize the projection matrix.  */


void myinit(void)
{
    
    /* attributes */
    
    glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */
    
    /* set up viewing */
    /* 512 x 512 window with origin lower left */
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 512.0, 0.0, 512.0);
    glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image. */

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */
    glRasterPos2i(0, 0);
    glDrawPixels(xmax_pixel, ymax_pixel, GL_RGB, GL_FLOAT, pixel_rgb);
    glFlush(); /* clear buffers */
}

/*  Main routine.  */

int main(int argc, char**argv)
{
    
    Read_Information();
    Setup_Parameters();
    Ray_Trace();
    
    /* Standard GLUT initialization */
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); /* default, not needed */
    glutInitWindowSize(1000, 1000); /* 500 x 500 pixel window */
    glutInitWindowPosition(0, 0); /* place window top left on display */
    glutCreateWindow("khindupur"); /* window title */
    glutDisplayFunc(display); /* display callback invoked when window opened */
    
    myinit(); /* set attributes */
    glutMainLoop(); /* enter event loop */
    
    return(0);
}