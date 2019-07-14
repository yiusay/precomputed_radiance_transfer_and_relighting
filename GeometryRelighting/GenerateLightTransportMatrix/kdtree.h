


struct KdNode;

typedef list<KdNode>::iterator KdNodeIter;

struct KdNode {

  union {
    int child_index;             // for inner node
    int triangles_header_index;  // for leaf node
  };

  int    split_dimension; // if split_dimension == -1 -> "leaf node"
  double split_coord;

  // for debug
  double bbmin[3], bbmax[3];
  // for debug
  
  KdNode(){
    split_dimension = -1;
  }
};


class Kdtree {

  double root_bbmin[3], root_bbmax[3]; // bbmin and bbmax of root node
  
  vector<KdNode> kd_nodes;
  vector< vector<int> > kdleaf_f_ids;

  void ConstructKdtree(int node_index, double *bbmin, double *bbmax,
                       vector<Vertex> &vertex, vector<Face> &face, vector<int> &kd_f_id, int level);
  void SetChildBBminBBmax(double *bbmin, double *bbmax, double (*child_bbmin)[3], double (*child_bbmax)[3],
                               double split_coord, int split_dimension);
  void DrawKdLeaf(double *bbmin, double *bbmax);

  bool ClipRaySegment(double *ray_origin, double *ray_direction, double *ray_direction_inv, double *bbmin, double *bbmax,
                      double &t_near, double &t_far);
  
public:
  Kdtree(){
    for(int i = 0; i < 3; i++){
      root_bbmin[i] =  1.0e6;
      root_bbmax[i] = -1.0e6;
    }
  }
  void InitKdtree(vector<Vertex> &vertex, vector<Face> &face);
  void DisplayKdtree();
  void DisplayLeafTriangles();

  bool GetRayIntersectingTriangles(double *ray_origin, double *ray_direction, int &closest_f_id, double& t_closest, double &u, double &v,
                                   vector<Vertex> &vertex, vector<Face> &face);
  void RenderUsingRayCasting(int width, int height, vector<Vertex> &vertex, vector<Face> &face);
  void RenderUsingPaintersAlgorithm(double *view_point, Mesh &face);
};

#define MAX_LEVEL 15

extern int intersect_triangle(double orig[3], double dir[3], double vert0[3],
                              double vert1[3], double vert2[3], double *t, double *u, double *v);


