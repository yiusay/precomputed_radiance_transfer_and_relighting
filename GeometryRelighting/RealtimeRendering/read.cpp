

#include "monteCarloPathTracing.h"
#include "mesh.h"


bool Mesh::ReadOFFFile(char *filename)
{
    FILE *fp;

    if((fp = fopen(filename, "r")) == NULL ){
        cerr << "file cannot be read.\n";
        return false;
    }

    char buf[512];

    // discard the first line
    fgets(buf, 512, fp);

    // obtain the number of vertices and faces
    fgets(buf, 512, fp);
    sscanf(buf, "%d%d", &n_vertices, &n_faces);

    int v_id = 0, f_id = 0;

    for(int i = 0; i < n_vertices; i++){
        double coord_in[3];

        fgets(buf, 512, fp);
        sscanf(buf, "%lf%lf%lf", &coord_in[0], &coord_in[1], &coord_in[2]);

        vertices.push_back( Vertex(coord_in, v_id++) );
    }

    for(int i = 0; i < n_faces; i++){
        int v_id[3], dummy;

        fgets(buf, 512, fp);
        sscanf(buf, "%d%d%d%d", &dummy, &v_id[0], &v_id[1], &v_id[2]);

        faces.push_back( Face(v_id, f_id++) );
    }

    cerr << "Reading Off file done...\n";

    fclose(fp);


    cerr << "# of vertices  " << n_vertices << endl;
    cerr << "# of faces     " << n_faces    << endl;


    double range_min[3] = {  1.0e6,  1.0e6,  1.0e6, };
    double range_max[3] = { -1.0e6, -1.0e6, -1.0e6, };
    double center[3];

     for(unsigned int i = 0; i < vertices.size(); i++){
        for(int j = 0; j < 3; j++){
            if(vertices[i].point[j] < range_min[j])	range_min[j] = vertices[i].point[j];
            if(vertices[i].point[j] > range_max[j])	range_max[j] = vertices[i].point[j];
        }
    }


    for(int i = 0; i < 3; i++) center[i] = (range_min[i] + range_max[i])*0.5;

    double largest_range = -1.0;

    for(int i = 0; i < 3; i++){
        if(largest_range < range_max[i]-range_min[i]) largest_range = range_max[i]-range_min[i];
    }

    double scale_factor = 400.0/largest_range;
    double min_y = 1.0e6;

    for(unsigned int i = 0; i < vertices.size(); i++){
        for(int j = 0; j < 3; j++){
            vertices[i].point[j] = (vertices[i].point[j] - center[j]) * scale_factor;
        }

        if(vertices[i].point[1] < min_y) min_y = vertices[i].point[1];
    }

    for(unsigned int i = 0; i < vertices.size(); i++){
        vertices[i].point[0] += 256.0;
        vertices[i].point[1] -= min_y;
        vertices[i].point[2] += 256.0;
    }

 
    for(unsigned int i = 0; i < faces.size(); i++) AssignFaceNormal(i);
    AssignVertexNormal();
    
    return true;
}




void Mesh::AssignFaceNormal(int f_id)
{
    double vec1[3], vec2[3];

  //      cerr << f_id << "  " << vertices.size() << " " << faces[f_id].v_id[0] << " " << faces[f_id].v_id[1] << " " << faces[f_id].v_id[2] << endl;
    for(int j = 0; j < 3; j++){
        vec1[j] = vertices[ faces[f_id].v_id[1] ].point[j] - vertices[ faces[f_id].v_id[0] ].point[j];
        vec2[j] = vertices[ faces[f_id].v_id[2] ].point[j] - vertices[ faces[f_id].v_id[0] ].point[j];
    }

    CrossProduct(vec1, vec2, faces[f_id].normal);
    GetArea(faces[f_id].normal, faces[f_id].area);
    Normalize(faces[f_id].normal);
}

void Mesh::AssignVertexNormal()
{
    // store face iteretors incident to each vertex
    vector< vector<int> > Ring(vertices.size());

    for(unsigned int i = 0; i < faces.size(); i++){
        Ring[faces[i].v_id[0]].push_back(i);
        Ring[faces[i].v_id[1]].push_back(i);
        Ring[faces[i].v_id[2]].push_back(i);
    }

    for(unsigned int i = 0; i < vertices.size(); i++){
        vertices[i].normal[0] = vertices[i].normal[1] = vertices[i].normal[2] = 0.0;
        vertices[i].area = 0.0;

        // traverse faces incident to "vi" in CCW
         for(unsigned int j = 0; j < Ring[i].size(); j++){

            vertices[i].normal[0] += faces[Ring[i][j]].normal[0] * faces[Ring[i][j]].area;
            vertices[i].normal[1] += faces[Ring[i][j]].normal[1] * faces[Ring[i][j]].area;
            vertices[i].normal[2] += faces[Ring[i][j]].normal[2] * faces[Ring[i][j]].area;
           
            vertices[i].area += faces[Ring[i][j]].area;
         }

        double invVertexArea = 1.0 / vertices[i].area;

        vertices[i].normal[0] *= invVertexArea;
        vertices[i].normal[1] *= invVertexArea;
        vertices[i].normal[2] *= invVertexArea;
    }

}


