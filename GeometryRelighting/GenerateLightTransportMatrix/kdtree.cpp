
#include "monteCarloPathTracing.h"
#include "mesh.h"
#include "kdtree.h"


//template<class T> inline void Swap(T &a, T &b)
//{
//  T t = a;  a = b;  b = t;
//}

template<class T>
inline void GetMinMaxAmongThreeValues(T v1, T v2, T v3, T &min, T &max)
{
  if(v1 < v2){
    if(v1 < v3){
      min = v1;
      if(v2 > v3) max = v2;
      else        max = v3;
      return;
    }
  }else{ // v2 < v1
    if(v2 < v3){
      min = v2;
      if(v1 > v3) max = v1;
      else        max = v3;
      return;
    }
  }

  min = v3;
  if(v1 > v2) max = v1;
  else        max = v2;
}

bool IsTriangleInsideCube(double *tri_point[3], double *tri_normal , double *bbmin, double *bbmax);
void GetProjectionCubeMinMax(double *projection_vector, double *bbmin, double *bbmax, double &min, double &max);


void Kdtree::InitKdtree(vector<Vertex> &vertex, vector<Face> &face)
{
  for(unsigned int i = 0; i < face.size(); i++){
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        if(vertex[ face[i].v_id[j] ].point[k] < root_bbmin[k]) root_bbmin[k] = vertex[ face[i].v_id[j] ].point[k];
        if(vertex[ face[i].v_id[j] ].point[k] > root_bbmax[k]) root_bbmax[k] = vertex[ face[i].v_id[j] ].point[k];
      }
    }
  }

  vector<int> kd_f_id(face.size());
  for(unsigned int i = 0; i < face.size(); i++) kd_f_id[i] = i;

  for(int i = 0; i < 3; i++){
    root_bbmin[i] -= (root_bbmax[i]-root_bbmin[i]) * 0.01;
    root_bbmax[i] += (root_bbmax[i]-root_bbmin[i]) * 0.01;
  }
  
  kd_nodes.push_back( KdNode() ); // as a root node
  
  ConstructKdtree(0, root_bbmin, root_bbmax, vertex, face, kd_f_id, 0);
}

void Kdtree::ConstructKdtree(int node_index, double *bbmin, double *bbmax,
                             vector<Vertex> &vertex, vector<Face> &face, vector<int> &kd_f_id, int level)
{
  if(level >= MAX_LEVEL){ // stop subdividing further

    kd_nodes[node_index].triangles_header_index = kdleaf_f_ids.size();
    kdleaf_f_ids.push_back(kd_f_id);
    
    
    // for debug
    for(int i = 0; i < 3; i++){
      kd_nodes[node_index].bbmin[i] = bbmin[i];
      kd_nodes[node_index].bbmax[i] = bbmax[i];
    }
    // for debug
    
    return;
    
  }else{

    double minimum_cost = 1.0e6, split_coord;
    int    split_dimension;

    double w = bbmax[0]-bbmin[0], h = bbmax[1]-bbmin[1], d = bbmax[2]-bbmin[2];
    double inv_surface_area = 1.0 / (w*h + h*d + d*w); // half of surface area of this node

    for(int i = 0; i < 3; i++){ // test spliting plane perpendicular to x, y, and z-axis, respectively

      vector<double> in_vertices;   // vertices of each triangle which has the smallest coordinate
      vector<double> out_vertices;  // vertices of each triangle which has the largest coordinate

      for(unsigned int j = 0; j < kd_f_id.size(); j++){
        Face *fp = &(face[ kd_f_id[j] ]);
        double smallest_coord, largest_coord;
        
        GetMinMaxAmongThreeValues(vertex[ fp->v_id[0] ].point[i],
                                  vertex[ fp->v_id[1] ].point[i],
                                  vertex[ fp->v_id[2] ].point[i],
                                  smallest_coord, largest_coord);

        in_vertices.push_back(smallest_coord);
        out_vertices.push_back(largest_coord);
      } // for(unsigned int j = 0; j < kd_f_id.size(); j++){

      sort(in_vertices.begin(), in_vertices.end());
      sort(out_vertices.begin(), out_vertices.end());

      unsigned int in_index = 0, out_index = 0;
      int num_child_faces[2] = { 0, kd_f_id.size() };

      while( in_index < in_vertices.size() || out_index < out_vertices.size() ){

        bool is_in_target;
        
        if(in_index >= in_vertices.size())                       is_in_target = false;
        else if(out_index >= out_vertices.size())                is_in_target = true;
        else if(in_vertices[in_index] < out_vertices[out_index]) is_in_target = true;  // choose the smallest one as a test split_coord
        else                                                     is_in_target = false;

        double target_coord;
        
        if(is_in_target){
          target_coord = in_vertices[in_index];
          in_index++;

          num_child_faces[0]++;
        }else{
          target_coord = out_vertices[out_index];
          out_index++;

          num_child_faces[1]--;
        }

        while( in_index < in_vertices.size() && fabs(target_coord - in_vertices[in_index]) < EPSILON ){
          num_child_faces[0]++;
          in_index++;
        }

        while( out_index < out_vertices.size() && fabs(target_coord - out_vertices[out_index]) < EPSILON ){
          num_child_faces[1]--;
          out_index++;
        }

        double surface_area[2];

        switch(i){
        case 0:
          surface_area[0] = (h*d + (target_coord-bbmin[0])*(h+d));
          surface_area[1] = (h*d + (bbmax[0]-target_coord)*(h+d));
          break;
        case 1:
          surface_area[0] = (d*w + (target_coord-bbmin[1])*(d+w));
          surface_area[1] = (d*w + (bbmax[1]-target_coord)*(d+w));
          break;
        case 2:
          surface_area[0] = (w*h + (target_coord-bbmin[2])*(w+h));
          surface_area[1] = (w*h + (bbmax[2]-target_coord)*(w+h));
          break;
        }

        double current_cost = (surface_area[0]*num_child_faces[0] + surface_area[1]*num_child_faces[1]) * inv_surface_area;

        if(current_cost < minimum_cost){
          minimum_cost = current_cost;
          split_coord = target_coord;
          split_dimension = i;
        }

      } // while( in_index < in_vertices.size() || out_index < out_vertices.size() ){

    } // for(int i = 0; i < 3; i++){


    if(minimum_cost > (double)kd_f_id.size()){ // we should not split further because the cost without splitting is smaller

      kd_nodes[node_index].triangles_header_index = kdleaf_f_ids.size();
      kdleaf_f_ids.push_back(kd_f_id);
      
      return;

    }else{

      vector<int> child_kd_f_id[2];
      double child_bbmin[2][3], child_bbmax[2][3];

      SetChildBBminBBmax(bbmin, bbmax, child_bbmin, child_bbmax, split_coord, split_dimension);

      for(unsigned int i = 0; i < kd_f_id.size(); i++){
        int *v_id = face[ kd_f_id[i] ].v_id;
        
        if( vertex[ v_id[0] ].point[split_dimension] < split_coord &&
            vertex[ v_id[1] ].point[split_dimension] < split_coord &&
            vertex[ v_id[2] ].point[split_dimension] < split_coord ){
          child_kd_f_id[0].push_back(kd_f_id[i]);

        }else if( vertex[ v_id[0] ].point[split_dimension] > split_coord &&
                  vertex[ v_id[1] ].point[split_dimension] > split_coord &&
                  vertex[ v_id[2] ].point[split_dimension] > split_coord ){
          child_kd_f_id[1].push_back(kd_f_id[i]);

        }else{
          double *triangle_point[3] = { vertex[ v_id[0] ].point, vertex[ v_id[1] ].point, vertex[ v_id[2] ].point };

          if( IsTriangleInsideCube(triangle_point, face[ kd_f_id[i] ].normal, child_bbmin[0], child_bbmax[0]) ){
            child_kd_f_id[0].push_back(kd_f_id[i]);
          }
          if( IsTriangleInsideCube(triangle_point, face[ kd_f_id[i] ].normal, child_bbmin[1], child_bbmax[1]) ){
            child_kd_f_id[1].push_back(kd_f_id[i]);
          }
        }

      } //  for(unsigned int i = 0; i < kd_f_id.size(); i++){

      int child_index[2] = { kd_nodes.size(),  kd_nodes.size()+1 };

      kd_nodes.push_back( KdNode() );
      kd_nodes.push_back( KdNode() );
        
      kd_nodes[node_index].child_index = child_index[0]; // we can access child_index[1] just by incrementing node->child_index, so no need to explicitely store it.
      kd_nodes[node_index].split_coord = split_coord;
      kd_nodes[node_index].split_dimension = split_dimension;

      ConstructKdtree(child_index[0], child_bbmin[0], child_bbmax[0], vertex, face, child_kd_f_id[0], level+1);
      ConstructKdtree(child_index[1], child_bbmin[1], child_bbmax[1], vertex, face, child_kd_f_id[1], level+1);
    } // else{

  } // else{

}



inline void Kdtree::SetChildBBminBBmax(double *bbmin, double *bbmax, double (*child_bbmin)[3], double (*child_bbmax)[3],
                                       double split_coord, int split_dimension)
{
  for(int i = 0; i < 3; i++){
    child_bbmin[0][i] = bbmin[i];
    child_bbmax[1][i] = bbmax[i];
  }
  
  switch(split_dimension){
  case 0:
    child_bbmax[0][0] = split_coord;
    child_bbmax[0][1] = bbmax[1];
    child_bbmax[0][2] = bbmax[2];
    child_bbmin[1][0] = split_coord;
    child_bbmin[1][1] = bbmin[1];
    child_bbmin[1][2] = bbmin[2];
    break;
  case 1:
    child_bbmax[0][0] = bbmax[0];
    child_bbmax[0][1] = split_coord;
    child_bbmax[0][2] = bbmax[2];
    child_bbmin[1][0] = bbmin[0];
    child_bbmin[1][1] = split_coord;
    child_bbmin[1][2] = bbmin[2];
    break;
  case 2:
    child_bbmax[0][0] = bbmax[0];
    child_bbmax[0][1] = bbmax[1];
    child_bbmax[0][2] = split_coord;
    child_bbmin[1][0] = bbmin[0];
    child_bbmin[1][1] = bbmin[1];
    child_bbmin[1][2] = split_coord;
    break;
  }

}



// intersection test between cube and triangle "Separating Axes"
bool IsTriangleInsideCube(double *tri_point[3], double *tri_normal , double *bbmin, double *bbmax)
{
  double vec[3][3];

  for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) vec[i][j] = tri_point[i][j] - bbmin[j];
  

  // "Separating Axes" using normal vector of triangle

  double dot = DotProduct(vec[0], tri_normal);

  double min, max;
  GetProjectionCubeMinMax(tri_normal, bbmin, bbmax, min, max);
  
  if(dot > max || dot < min) return false;


  // "Separating Axes" using 3 normal vectors of cube
  for(int i = 0; i < 3; i++){

    // when i == 0, projection to quad normal = (1, 0, 0)
    // when i == 1, projection to quad normal = (0, 1, 0)
    // when i == 2, projection to quad normal = (0, 0, 1)

    double dot[3];
    for(int j = 0; j < 3; j++) dot[j] = vec[j][i];
    
    double dot_min, dot_max;
    GetMinMaxAmongThreeValues(dot[0], dot[1], dot[2], dot_min, dot_max);
    
    if(dot_min > bbmax[i]-bbmin[i] || dot_max < 0.0) return false;
  }


  // "Separating Axes" using 9 pairs of cross products
  double TriVec[3], CrossVec[3];

  for(int i = 0; i < 3; i++){
    
    for(int j = 0; j < 3; j++) TriVec[j] = tri_point[(i+1)%3][j] - tri_point[i][j];

    for(int j = 0; j < 3; j++){

      switch(j){
      case 0:
        CrossVec[0] = 0.0;
        CrossVec[1] = -TriVec[2];
        CrossVec[2] = TriVec[1]; 
        break;
      case 1:
        CrossVec[0] = TriVec[2];
        CrossVec[1] = 0.0;
        CrossVec[2] = -TriVec[0]; 
        break;
      case 2:
        CrossVec[0] = -TriVec[1];
        CrossVec[1] = TriVec[0];
        CrossVec[2] = 0.0; 
        break;
      }

      double dot[3], dot_min, dot_max;

      for(int j = 0; j < 3; j++) dot[j] = DotProduct(vec[j], CrossVec);

      GetMinMaxAmongThreeValues(dot[0], dot[1], dot[2], dot_min, dot_max);
      
      double min, max;
      GetProjectionCubeMinMax(CrossVec, bbmin, bbmax, min, max);

      if(dot_min > max || dot_max < min) return false;
    } // for(int j = 0; j < 3; j++){

  } // for(int i = 0; i < 3; i++){

  
  // The given triangle is intersecting with the cube
  return true;
}



void GetProjectionCubeMinMax(double *projection_vector, double *bbmin, double *bbmax, double &min, double &max)
{
  // min = 0.0 at maximum, max = 0.0 at minimum
  min = max = 0.0;

  for(int i = 0; i <= 6; i++){
    double vec[3];
    
    switch(i){
    case 0:  // vector fomr bbmin to vertex with id "0"isame for belowj
      vec[0] = bbmax[0] - bbmin[0];
      vec[1] = bbmax[1] - bbmin[1];
      vec[2] = bbmax[2] - bbmin[2];
      break;
    case 1:  // "1"
      vec[0] = 0.0;
      vec[1] = bbmax[1] - bbmin[1];
      vec[2] = bbmax[2] - bbmin[2];
      break;
    case 2:  // "2"
      vec[0] = bbmax[0] - bbmin[0];
      vec[1] = 0.0;
      vec[2] = bbmax[2] - bbmin[2];
      break;
    case 3:  // "3"
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = bbmax[2] - bbmin[2];
      break;
    case 4:  // "4"
      vec[0] = bbmax[0] - bbmin[0];
      vec[1] = bbmax[1] - bbmin[1];
      vec[2] = 0.0;
      break;
    case 5:  // "5"
      vec[0] = 0.0;
      vec[1] = bbmax[1] - bbmin[1];
      vec[2] = 0.0;
      break;
    case 6:  // "6"
      vec[0] = bbmax[0] - bbmin[0];
      vec[1] = 0.0;
      vec[2] = 0.0;
      break;
    }

    double dot = DotProduct(projection_vector, vec);
    if(dot > max) max = dot;
    if(dot < min) min = dot;
    
  } // for(int i = 0; i <= 6; i++){
  
}



struct NodeData {
  int    kdnode_index;
  double bbmin[3], bbmax[3];

  NodeData(int in_kdnode_index, double *in_bbmin, double *in_bbmax){
    kdnode_index = in_kdnode_index;
    for(int i = 0; i < 3; i++){
      bbmin[i] = in_bbmin[i];
      bbmax[i] = in_bbmax[i];
    }
  }
};



bool Kdtree::ClipRaySegment(double *ray_origin, double *ray_direction, double *ray_direction_inv, double *bbmin, double *bbmax,
                            double &t_near, double &t_far)
{
  t_near = EPSILON;
  t_far  = 1.0e6;

  for(int i = 0; i < 3; i++){
    
    if( fabs(ray_direction_inv[i]) > EPSILON ){
      
      double t1 = (bbmin[i]-ray_origin[i]) * ray_direction_inv[i];
      double t2 = (bbmax[i]-ray_origin[i]) * ray_direction_inv[i];      
      
      if(ray_direction_inv[i] > 0.0){
        if(t1 > t_near) t_near = t1;
        if(t2 < t_far)  t_far  = t2;
      }else{
        if(t1 < t_far)  t_far  = t1;        
        if(t2 > t_near) t_near = t2;
      }
      
    } // if( fabs(ray_direction_inv[i]) > EPSILON ){

  } // for(int i = 0; i < 3; i++){

	return true;
}


struct NodeAndRayData {
  int    node_id;
  double t_near, t_far;

  NodeAndRayData(int node_id_in, double t_near_in, double t_far_in){
    node_id = node_id_in;
    t_near  = t_near_in;
    t_far   = t_far_in;
  }
};







bool Kdtree::GetRayIntersectingTriangles(double *ray_origin, double *ray_direction, int &closest_f_id, 
                                         double& t_closest, double &u_closest, double &v_closest, 
                                         vector<Vertex> &vertex, vector<Face> &face)
{
  double t_near, t_far;
  double ray_direction_inv[3];

  for(int i = 0; i < 3; i++){
    if( fabs(ray_direction[i]) > EPSILON ) ray_direction_inv[i] = 1.0/ray_direction[i];
    else                                   ray_direction_inv[i] = 0.0;
  }

  if( ClipRaySegment(ray_origin, ray_direction, ray_direction_inv, root_bbmin, root_bbmax, t_near, t_far) == false) return false;

  if(t_near > t_far) return false; // ray does not intersect with root node 
  
  
  stack<NodeAndRayData> node_ray_data;

  int current_node_id = 0;

  closest_f_id = -1;
  t_closest    = 1.0e6;
  
  while(true){

    while(kd_nodes[current_node_id].split_dimension != -1){ // if split_dimension == -1 -> "leaf node"

      int    dim = kd_nodes[current_node_id].split_dimension;
      double d = (kd_nodes[current_node_id].split_coord - ray_origin[dim]) * ray_direction_inv[dim];
      
      if( fabs(ray_direction_inv[dim]) < EPSILON ){ // ray is parallel to the splitting plane
        
        if(ray_origin[dim] < kd_nodes[current_node_id].split_coord)
          current_node_id = kd_nodes[current_node_id].child_index;
        else
          current_node_id = kd_nodes[current_node_id].child_index + 1; 
        
      }else if (d <= t_near) {
        // case one, d <= t_near <= t_far -> cull front side

        if(ray_direction[dim] < 0.0) current_node_id = kd_nodes[current_node_id].child_index;
        else                         current_node_id = kd_nodes[current_node_id].child_index + 1; 

      } else if (d >= t_far) {
        // case two, t_near <= t_far <= d -> cull back side

        if(ray_direction[dim] > 0.0) current_node_id = kd_nodes[current_node_id].child_index;
        else                         current_node_id = kd_nodes[current_node_id].child_index + 1;
        
      } else {
        // case three: traverse both sides in turn

        if(ray_direction[dim] > 0.0){
          node_ray_data.push( NodeAndRayData(kd_nodes[current_node_id].child_index+1, d, t_far) );
          current_node_id = kd_nodes[current_node_id].child_index;
          t_far = d;
        }else{
          node_ray_data.push( NodeAndRayData(kd_nodes[current_node_id].child_index, d, t_far) );
          current_node_id = kd_nodes[current_node_id].child_index+1;
          t_far = d;
        }

      }
      
    } // while(kd_nodes[current_node_id].split_dimension != -1){

    for(unsigned int i = 0; i < kdleaf_f_ids[ kd_nodes[current_node_id].triangles_header_index ].size(); i++){

      int f_id = kdleaf_f_ids[ kd_nodes[current_node_id].triangles_header_index ][i];
      
      double *triangle_point[3] = { vertex[ face[f_id].v_id[0] ].point, vertex[ face[f_id].v_id[1] ].point, vertex[ face[f_id].v_id[2] ].point };
      
      double t, u, v;
      if( intersect_triangle(ray_origin, ray_direction,
                             triangle_point[0], triangle_point[1], triangle_point[2],
                             &t, &u, &v) ){

        if(t > EPSILON){
          if(t < t_closest){
            closest_f_id = f_id;
            t_closest = t;
          }
        }
          
          
      } // if( intersect_triangle(...) ){
      
    } // for(unsigned int i = 0; i < kdleaf_triangles[ kd_nodes[current_node_id].triangles_header_index ].size(); i++){

    // nothing else to traverse any more
    if(node_ray_data.empty()) break;
    
    current_node_id = ( node_ray_data.top() ).node_id;
    t_near          = ( node_ray_data.top() ).t_near;
    t_far           = ( node_ray_data.top() ).t_far;

    node_ray_data.pop();
    
  } // while(true){
  
  
  if(closest_f_id != -1) return true;
  else                   return false;
}

