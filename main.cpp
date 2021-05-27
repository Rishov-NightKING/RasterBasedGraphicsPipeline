#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <stack>
#include <random>
#include <limits>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;

class Point{

public:
    double x, y, z, w;

    Point()
    {
        this->x = this->y = this->z = 0.0;
        this->w = 1;
    }
    Point(double x, double y, double z) : x(x), y(y), z(z), w(1) {}

    void normalize_point()
    {
        double value = sqrt(x * x + y * y + z * z);
        x /= value;
        y /= value;
        z /= value;
    }

    void printPoint() const
    {
        cout << "(" << setprecision(7) << fixed << this->x << ", "
            <<setprecision(7) << fixed << this->y << ", "
            << setprecision(7) << fixed << this->z << ")" << endl;
    }
};

class Matrix{
    int row, column;
public:
    vector< vector<double> > mat;
    Matrix() {
        //default 4X4 matrix
        this->row = this->column = 4;
        initialize_matrix();
    }

    Matrix(int row, int column) {
        this->row = row;
        this->column = column;
        initialize_matrix();
    }


    int getRow() const {
        return row;
    }

    int getColumn() const {
        return column;
    }

    void initialize_matrix()
    {
        mat.resize(row);
        for(int i = 0; i < row; i++)
        {
            mat[i].resize(column);
            for(int j = 0; j < column; j++)
            {
                mat[i][j] = 0.0;
            }
        }
    }

    void convert_to_Identity_matrix()
    {
        if(this->row == this->column)
        {
            for(int i = 0; i < row; i++)
            {
                for(int j = 0; j < column; j++)
                {
                    if(i == j) mat[i][j] = 1.0;
                    else mat[i][j] = 0.0;
                }
            }
        }
        else cout << "Cannot convert to identity matrix as rows and columns are not equal" << endl;
    }

    void scale_matrix_by_column_last()
    {
        for(int i = 0; i < column - 1; i++)
        {
            double column_last = mat[row - 1][i];
            if(column_last != 1)
            {
                for(int j = 0; j < row; j++)
                {
                    mat[j][i] /= column_last;
                }
            }
        }
    }

    void print_matrix()
    {
        for(int i = 0; i < row; i++)
        {
            for(int j = 0; j < column; j++)
            {
                cout << setprecision(7) << fixed << mat[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_matrix_in_file(FILE *file)
    {
        //(n - 1) X (n - 1) vector matrix print in file
        for(int i = 0; i < row - 1; i++)
        {
            for(int j = 0; j < column - 1; j++)
            {
                fprintf(file, "%.7lf", mat[j][i]);
                if(j < column - 2) fprintf(file, " ");
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n");
    }

};

//C++ Modern approach
int random_in_integer_range(int range_from, int range_to)
{
    std::random_device rand_dev; // obtain a random number from hardware
    std::mt19937 generator(rand_dev()); // seed the generator
    std::uniform_int_distribution<> distribution(range_from, range_to); // define the range.. it is inclusive
    return distribution(generator);
}

double random_in_double_range(double range_from, double range_to)
{
    std::random_device rand_dev; // obtain a random number from hardware
    std::mt19937 generator(rand_dev()); // seed the generator
    std::uniform_real_distribution<> distribution(range_from, range_to); // define the range.. it is inclusive
    return distribution(generator);
}

class Triangle{

public:
    vector<Point> end_points;
    vector<int> RGB_color;
    double max_X, min_X, max_Y, min_Y;

    Triangle()
    {
        end_points.resize(3);
        RGB_color.resize(3);
        for(int i = 0; i < 3; i++)
        {
            end_points[i] = Point();
            RGB_color[i] = 0;
        }
        max_X = max_Y = numeric_limits<double>::min();
        min_X = min_Y = numeric_limits<double>::max();
    }


    void set_random_color()
    {
        for(int i = 0; i < 3; i++)
        {
            RGB_color[i] = random_in_integer_range(0, 255);
        }
    }

    void print_triangle()
    {
        cout << "\nTriangle Info" << endl;
        for(int i = 0; i < 3; i++)
        {
            cout << "Endpoint-" << (i + 1) << ":  ";
            end_points[i].printPoint();
        }
        cout << "RGB color value:  R: " << RGB_color[0] << "   G: " << RGB_color[1] << "   B: " << RGB_color[2] << endl;
        cout << "min_X: " << min_X << "  max_X: " << max_X << endl;
        cout << "min_Y: " << min_Y << "  max_Y: " << max_Y << endl;
    }

};
/************************GLOBAL VARIABLES***************************************/
int triangles = 0;

Point eye, look, up;
double fovY, aspectRatio, near, far;

int Screen_Width, Screen_Height;
double left_limit_X, right_limit_X;
double bottom_limit_Y, top_limit_Y;
double front_limit_Z, rear_limit_Z;

stack<Matrix> Stack;
Matrix V, P;

/******************************************************************************/



double degreeToRadianAngle(double degree)
{
    return pi / 180 * degree;
}
double vector_dot_product(Point a, Point b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Point vector_cross_product(Point a, Point b)
{
    Point temp;
    temp.x = a.y * b.z - a.z * b.y;
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;
    return temp;
}
Matrix matrix_multiply(Matrix a, Matrix b)
{
    if(a.getColumn() != b.getRow())
    {
        cout<<"Matrix multiplication cannot be performed."<<endl;
    }
    else
    {
        int mul_row = a.getRow();
        int mul_col = b.getColumn();
        Matrix c(mul_row, mul_col);
        for(int i=0;i<mul_row;i++){
            for(int j=0;j<mul_col;j++){
                double sum = 0;
                for(int k=0;k<a.getColumn();k++){
                    sum += a.mat[i][k] * b.mat[k][j];
                }
                c.mat[i][j] = sum;
            }
        }
        return c;
    }
    return Matrix();
}
void operation_on_matrix(Matrix mat)
{
    if(Stack.empty()) cout << "Error during modeling transformation....Stack empty" <<endl;
    else
    {
        Matrix stack_top = matrix_multiply(Stack.top(), mat);
        Stack.pop();
        Stack.push(stack_top);
    }
}
Point Rodrigues_rotation(Point x, Point a, double angle)
{
    double sin_theta = sin(degreeToRadianAngle(angle));
    double cos_theta = cos(degreeToRadianAngle(angle));
    double temp = (1 - cos_theta) * vector_dot_product(a, x);

    Point C = vector_cross_product(a, x);

    Point R;
    R.x = cos_theta * x.x + temp * a.x + sin_theta * C.x;
    R.y = cos_theta * x.y + temp * a.y + sin_theta * C.y;
    R.z = cos_theta * x.z + temp * a.z + sin_theta * C.z;

    return R;
}
void gluLookAt()
{
    cin >> eye.x >> eye.y >> eye.z;
    cin >> look.x >> look.y >> look.z;
    cin >> up.x >> up.y >> up.z;
}
void gluPerspective()
{
    cin >> fovY >> aspectRatio >> near >> far;
}
void view_transformation()
{
    //View Transformation
    Point l = Point(look.x - eye.x, look.y - eye.y, look.z - eye.z);
    l.normalize_point();

    Point r = vector_cross_product(l, up);
    r.normalize_point();

    Point u = vector_cross_product(r, l);

    //T
    Matrix T;
    T.convert_to_Identity_matrix();
    T.mat[0][3] = -eye.x;
    T.mat[1][3] = -eye.y;
    T.mat[2][3] = -eye.z;

    Matrix R;
    R.mat[0][0] = r.x;
    R.mat[0][1] = r.y;
    R.mat[0][2] = r.z;

    R.mat[1][0] = u.x;
    R.mat[1][1] = u.y;
    R.mat[1][2] = u.z;

    R.mat[2][0] = -l.x;
    R.mat[2][1] = -l.y;
    R.mat[2][2] = -l.z;

    R.mat[3][3] = 1;

    V = matrix_multiply(R, T);
}
void projection_transformation()
{
    //Projection Transformation
    double fovX = fovY * aspectRatio;
    double t_projection = near * tan(degreeToRadianAngle(fovY / 2));
    double r_projection = near * tan(degreeToRadianAngle(fovX / 2));

    P.mat[0][0] = near / r_projection;
    P.mat[1][1] = near / t_projection;
    P.mat[2][2] = -(far + near) / (far - near);
    P.mat[2][3] = -(2 * far * near) / (far - near);
    P.mat[3][2] = -1;
}

double distance_between_points(Point a, Point b)
{
    return sqrt( pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

vector<Point> determine_intersecting_points(Point A, Point B, Point C, Triangle triangle)
{
    vector<Point> pair_vector;
    
    if(distance_between_points(triangle.end_points[0], A) + distance_between_points(triangle.end_points[1], A) - distance_between_points(triangle.end_points[0], triangle.end_points[1]) < 0.0001)
    {
        pair_vector.push_back(A);
    }
    if(distance_between_points(triangle.end_points[0], B) + distance_between_points(triangle.end_points[2], B) - distance_between_points(triangle.end_points[0], triangle.end_points[2]) < 0.0001)
    {
        pair_vector.push_back(B);
    }
    if(distance_between_points(triangle.end_points[1], C) + distance_between_points(triangle.end_points[2], C) - distance_between_points(triangle.end_points[1], triangle.end_points[2]) < 0.0001)
    {
        pair_vector.push_back(C);
    }
    return pair_vector;
}
void print_vector_in_file(FILE *file, vector< vector<double> > vec)
{
    for(int i = 0; i < Screen_Height ; i++)
    {
        for(int j = 0; j < Screen_Width ; j++)
        {
            if(vec[i][j] < rear_limit_Z) fprintf(file, "%.6lf", vec[i][j]);
            if(j < Screen_Width - 1) fprintf(file, "  ");
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}
int main() {
    Matrix matrix, triangle_vector, mat_proj, mat_view;
    double x, y, z, angle;
    double dx, dy;
    double top_Y, bottom_Y, left_X, right_X;
    string command;
    vector<Triangle> triangle_store;

    /****************************OPEN FILES***************************************/
    FILE *scene = freopen("scene.txt", "r", stdin);
    FILE *stage1 = fopen("stage1.txt","w");
    FILE *stage2 = fopen("stage2.txt","w");
    FILE *stage3 = fopen("stage3.txt","w");
    FILE *z_buffer_file = fopen("z_buffer_1605084.txt", "w");
    /****************************************************************************/


    gluLookAt();
    gluPerspective();
    view_transformation();
    projection_transformation();

    //Identity matrix push
    matrix.convert_to_Identity_matrix();
    Stack.push(matrix);

    while (true)
    {
        cin >> command;
        if(command == "end") break;
        else if(command == "triangle")
        {
            triangles++;
            for(int i = 0; i < triangle_vector.getRow() - 1; i++){
                for(int j = 0; j < triangle_vector.getColumn() - 1; j++)
                {
                    cin >> triangle_vector.mat[j][i];
                }
            }
            for(int i = 0; i < triangle_vector.getRow() - 1; i++)
            {
                triangle_vector.mat[3][i] = 1;
                //triangle_vector.mat[i][3] = 1;
            }

            if(Stack.empty())
            {
                cout << "Triangle operation error" <<endl;
            }
            else
            {
                triangle_vector = matrix_multiply(Stack.top(), triangle_vector);
                triangle_vector.print_matrix_in_file(stage1);

                //View transformation
                mat_view = matrix_multiply(V, triangle_vector);
                mat_view.print_matrix_in_file(stage2);

                //Projection Transformation
                mat_proj = matrix_multiply(P, mat_view);
                mat_proj.scale_matrix_by_column_last();
                mat_proj.print_matrix_in_file(stage3);
            }

        }
        else if(command == "translate")
        {
            matrix.convert_to_Identity_matrix();
            for(int i=0;i<3;i++)
            {
                cin >> matrix.mat[i][3];
            }

            operation_on_matrix(matrix);

        }
        else if(command == "scale")
        {
            matrix.convert_to_Identity_matrix();
            for(int i=0;i<3;i++)
            {
                cin >> matrix.mat[i][i];
            }

            operation_on_matrix(matrix);
        }
        else if(command == "rotate")
        {
            cin >> angle >> x >> y >> z;
            matrix.convert_to_Identity_matrix();

            //vector a = (x, y, z) and normalize a
            double value = sqrt(x * x + y * y + z * z);
            x /= value;
            y /= value;
            z /= value;

            Point i(1, 0, 0), j(0, 1, 0), k(0, 0, 1), a(x, y, z);

            Point C1 = Rodrigues_rotation(i, a, angle);
            Point C2 = Rodrigues_rotation(j, a, angle);
            Point C3 = Rodrigues_rotation(k, a, angle);

            matrix.mat[0][0] = C1.x;
            matrix.mat[1][0] = C1.y;
            matrix.mat[2][0] = C1.z;

            matrix.mat[0][1] = C2.x;
            matrix.mat[1][1] = C2.y;
            matrix.mat[2][1] = C2.z;

            matrix.mat[0][2] = C3.x;
            matrix.mat[1][2] = C3.y;
            matrix.mat[2][2] = C3.z;

            operation_on_matrix(matrix);
        }
        else if(command == "push")
        {
            if(!Stack.empty()) Stack.push(Stack.top());

        }
        else if(command == "pop")
        {
            if(!Stack.empty()) Stack.pop();
            else cout << "Pop not possible....Stack empty." <<endl;
        }
    }

    //READ config.txt file
    FILE *config = fopen("config.txt", "r");
    fscanf(config, "%d %d", &Screen_Width, &Screen_Height);
    fscanf(config, "%lf", &left_limit_X);
    fscanf(config, "%lf", &bottom_limit_Y);
    fscanf(config, "%lf %lf", &front_limit_Z, &rear_limit_Z);

    right_limit_X = left_limit_X * (-1);
    top_limit_Y = bottom_limit_Y * (-1);

    fclose(stage3);
    stage3 = fopen("stage3.txt","r");


    while (triangles--)
    {
        Triangle triangle;
        for(int i = 0; i < 3; i++)
        {
            fscanf(stage3, "%lf %lf %lf", &triangle.end_points[i].x, &triangle.end_points[i].y, &triangle.end_points[i].z);
            triangle.max_X = max(triangle.end_points[i].x, triangle.max_X);
            triangle.max_Y = max(triangle.end_points[i].y, triangle.max_Y);
            triangle.min_X = max(triangle.end_points[i].x, triangle.min_X);
            triangle.min_Y = max(triangle.end_points[i].y, triangle.min_Y);
        }
        triangle.set_random_color();
        triangle_store.push_back(triangle);
    }

    /* a.Create a pixel mapping between the x-y range values and the Screen_Width X
    Screen_height range.Besides, specify Top_Y and Left_X values. */
    dx = (right_limit_X - left_limit_X) / Screen_Width;
    dy = (top_limit_Y - bottom_limit_Y) / Screen_Height;
    top_Y = top_limit_Y - dy/2;
    bottom_Y = bottom_limit_Y + dy/2;
    left_X = left_limit_X + dx/2;
    right_X = right_limit_X - dx/2;


    /* c. Create a z-buffer, a two dimensional array of Screen_Width X Screen_Height dimension.
    Initialize all the values in z-buffer with z_max(far distance of eye) */
    vector< vector<double> > z_buffer(Screen_Height, vector<double>(Screen_Width, rear_limit_Z));


    /* d. Create a bitmap_image object with Screen_Width X Screen_Height resolution and initialize
    its background color with black.*/
    bitmap_image image(Screen_Width, Screen_Height);
    for(int i = 0; i < Screen_Height; i++)
    {
        for(int j = 0; j < Screen_Width; j++)
        {
            image.set_pixel(i , j, 0, 0, 0);
        }
    }


    int top_scanline_row, bottom_scanline_row, left_intersecting_col, right_intersecting_col;
    double clipping_min_Y, clipping_max_Y;

    for(int i = 0; i < triangle_store.size(); i++)
    {
       /* Find top_scanline and bottom_scanline after necessary clipping */
       clipping_max_Y = max(triangle_store[i].max_Y, top_Y);
       clipping_min_Y = min(triangle_store[i].min_Y, bottom_Y);

       /*  b. During scanning from top to bottom and left to right, check for the middle values of each
      cell. e.g. Top_Y- row_no*dy, Left_X + col_no*dx    */

       top_scanline_row = round((clipping_max_Y - bottom_Y) / dy);
       bottom_scanline_row = round((clipping_min_Y - bottom_Y) / dy);

       for(int ys = bottom_scanline_row; ys <= top_scanline_row; ys++)
       {
           double const1, const2, const3;
           double x1, x2, x3, y1, y2, y3, z1, z2, z3, xa, xb, xc, za, zb, zc, Za, Zb, Xa, Xb, zp, xp;
           x1 = triangle_store[i].end_points[0].x;
           x2 = triangle_store[i].end_points[1].x;
           x3 = triangle_store[i].end_points[2].x;

           y1 = triangle_store[i].end_points[0].y;
           y2 = triangle_store[i].end_points[1].y;
           y3 = triangle_store[i].end_points[2].y;

           z1 = triangle_store[i].end_points[0].z;
           z2 = triangle_store[i].end_points[1].z;
           z3 = triangle_store[i].end_points[2].z;


           double ysx = bottom_Y + ys * dy;
           const1 = (ysx - y1) / (y2 - y1);
           const2 = (ysx - y1) / (y3 - y1);
           const3 = (ysx - y2) / (y3 - y2);

           Point A, B, C;
           if(!isnan(const1) && !isinf(const1))
           {
               xa = x1 + const1 * (x2 - x1);
               za = z1 + const1 * (z2 - z1);
               A.x = xa;
               A.y = ysx;
               A.z = za;
           }

           if(!isnan(const2) && !isinf(const2))
           {
               xb = x1 + const2 * (x3 - x1);
               zb = z1 + const2 * (z3 - z1);
               B.x = xb;
               B.y = ysx;
               B.z = zb;
           }

           if(!isnan(const3) && !isinf(const3))
           {
               xc = x2 + const3 * (x3 - x2);
               zc = z2 + const3 * (z3 - z2);
               C.x = xc;
               C.y = ysx;
               C.z = zc;
           }

           vector<Point> intersecting_points = determine_intersecting_points(A, B, C, triangle_store[i]);

           if(intersecting_points.size() != 2) continue;

           Xa = intersecting_points[0].x;
           Xb = intersecting_points[1].x;
           Za = intersecting_points[0].z;
           Zb = intersecting_points[1].z;

           left_intersecting_col = round( (max(min(Xa, Xb), left_X ) - left_X) / dx);
           right_intersecting_col = round( (min(max(Xa, Xb), right_X)  - left_X) / dx);

           double const_zp = (Zb - Za) / (Xb - Xa);
           for (int xs = left_intersecting_col; xs <= right_intersecting_col; xs++)
           {
               xp = left_X + xs * dx;
               zp = Za + const_zp * (xp - Xa);

               if (zp < z_buffer[ys][xs] && zp > front_limit_Z)
               {
                   z_buffer[ys][xs] = zp;
                   image.set_pixel(xs, Screen_Height - ys, triangle_store[i].RGB_color[0], triangle_store[i].RGB_color[1],
                                   triangle_store[i].RGB_color[2]);
               }

           }
       }
    }
    print_vector_in_file(z_buffer_file, z_buffer);
    image.save_image("output_1605084.bmp");

    /****************************CLOSE FILES***************************************/

    fclose(config);
    fclose(scene);
    fclose(stage1);
    fclose(stage2);
    fclose(stage3);
    fclose(z_buffer_file);

    /****************************************************************************/
    return 0;
}
