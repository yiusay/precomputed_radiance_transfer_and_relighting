

struct Vertex {
    double point[3];
    double normal[3];
    double area;
    int id;

    Vertex(){};
    Vertex(double *data, int n){
        for(int i=0; i<3; i++) point[i] = data[i];
        id = n;
    }
    Vertex(double d_in0, double d_in1, double d_in2, int n){
        point[0] = d_in0;
        point[1] = d_in1;
        point[2] = d_in2;
        id = n;
    }
};


struct Face {
    int v_id[3];
    double centroid[3];

    double normal[3];
    double area;
    int id;

    Face(){}
    Face(int *data, int a){
        for(int i=0; i<3; i++) v_id[i] = data[i];
        id = a;
    }
    Face(int d_in0, int d_in1, int d_in2, int n){
        v_id[0] = d_in0;
        v_id[1] = d_in1;
        v_id[2] = d_in2;
        id = n;
    }
};


class Mesh{
protected:

    int n_vertices, n_faces;

    void AssignFaceNormal(int f_id);
    void AssignVertexNormal();

public:
    vector<Vertex> vertices;
    vector<Face>   faces;

    Mesh(){
        n_vertices = n_faces = 0;
    }

    bool ReadOFFFile(char *filename);
};
